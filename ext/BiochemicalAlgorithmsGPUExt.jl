module BiochemicalAlgorithmsGPUExt

# ============================================================================
#  GPU backend for the CompiledForceField (device-agnostic).
#
#  Loaded automatically when KernelAbstractions + Atomix are available (which
#  the Metal / CUDA package extensions pull in). A single kernel source thus
#  runs on any KernelAbstractions device; the Metal/CUDA extensions only need to
#  register their device via `_register_gpu_device!`.
#
#  Layout: scalar struct-of-arrays on the device (rx/ry/rz, Fx/Fy/Fz, pair
#  arrays). Nonbonded forces use atomic scatter-add; the cheap, irregular bonded
#  terms run on the host. The per-pair math is exactly the `_f_*_factor` /
#  `_e_*_pair` helpers used by the CPU loops, so results match within Float32
#  rounding.
# ============================================================================

using BiochemicalAlgorithms
using KernelAbstractions
using KernelAbstractions: @kernel, @index, @Const
import Atomix

import BiochemicalAlgorithms:
    CompiledForceField, EvalBackend, ForceField, Vector3,
    _GPU_DEVICE, _make_gpu_backend, _backend_init!, _backend_upload_pairs!,
    compute_energy!, compute_forces!,
    _finish_energy!, _force_bonded!, _apply_constraints!,
    _f_lj_factor, _f_hb_factor, _f_es_factor, _e_lj_pair, _e_hb_pair, _e_es_pair

const _KA = KernelAbstractions

struct KAGPUBackend{Dev} <: EvalBackend
    device::Dev
    workgroupsize::Int
end

function _make_gpu_backend(::ForceField, ::Type{Acc}) where Acc
    reg = _GPU_DEVICE[]
    reg === nothing && throw(ArgumentError(
        "the :gpu backend needs a GPU package loaded — run `using Metal` " *
        "(Apple Silicon) or `using CUDA`"))
    (Acc === Float64 && !reg.f64) && throw(ArgumentError(
        "this GPU device does not support Float64; use accumulation=Float32"))
    KAGPUBackend(reg.device, 256)
end

# device-resident mirror + host staging buffers (allocation-free per eval)
mutable struct GPUState{FV, IV, HV}
    rx::FV; ry::FV; rz::FV
    Fx::FV; Fy::FV; Fz::FV
    lj_i::IV; lj_j::IV; lj_A::FV; lj_B::FV; lj_s::FV; lj_e::FV
    hb_i::IV; hb_j::IV; hb_A::FV; hb_B::FV; hb_s::FV; hb_e::FV
    es_i::IV; es_j::IV; es_q::FV; es_s::FV; es_e::FV
    hx::HV; hy::HV; hz::HV         # host staging for positions
    hFx::HV; hFy::HV; hFz::HV      # host staging for forces
end

@inline _dev(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc} = cff.backend.device
@inline _wg(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc} = cff.backend.workgroupsize

_dalloc(dev, ::Type{E}, n) where E = _KA.allocate(dev, E, n)

# ---- (re)allocate + upload the nonbonded pair arrays -----------------------
function _upload_pairs!(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc}
    dev = _dev(cff)
    g = cff.gpu::GPUState
    up_i(host) = (d = _dalloc(dev, Int32, length(host)); copyto!(d, host); d)
    up_f(host) = (d = _dalloc(dev, Acc, length(host)); copyto!(d, Acc.(host)); d)
    g.lj_i = up_i(cff.lj.i); g.lj_j = up_i(cff.lj.j)
    g.lj_A = up_f(cff.lj.A); g.lj_B = up_f(cff.lj.B); g.lj_s = up_f(cff.lj.scaling)
    g.lj_e = _dalloc(dev, Acc, length(cff.lj.i))
    g.hb_i = up_i(cff.hb.i); g.hb_j = up_i(cff.hb.j)
    g.hb_A = up_f(cff.hb.A); g.hb_B = up_f(cff.hb.B); g.hb_s = up_f(cff.hb.scaling)
    g.hb_e = _dalloc(dev, Acc, length(cff.hb.i))
    g.es_i = up_i(cff.es.i); g.es_j = up_i(cff.es.j)
    g.es_q = up_f(cff.es.q1q2); g.es_s = up_f(cff.es.scaling)
    g.es_e = _dalloc(dev, Acc, length(cff.es.i))
    nothing
end

function _backend_init!(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc}
    dev = _dev(cff); n = cff.natoms
    FV = typeof(_dalloc(dev, Acc, 0)); IV = typeof(_dalloc(dev, Int32, 0))
    pos() = _dalloc(dev, Acc, n)
    cff.gpu = GPUState{FV, IV, Vector{Acc}}(
        pos(), pos(), pos(), pos(), pos(), pos(),
        _dalloc(dev, Int32, 0), _dalloc(dev, Int32, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0),
        _dalloc(dev, Int32, 0), _dalloc(dev, Int32, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0),
        _dalloc(dev, Int32, 0), _dalloc(dev, Int32, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0), _dalloc(dev, Acc, 0),
        Vector{Acc}(undef, n), Vector{Acc}(undef, n), Vector{Acc}(undef, n),
        Vector{Acc}(undef, n), Vector{Acc}(undef, n), Vector{Acc}(undef, n),
    )
    _upload_pairs!(cff)
    nothing
end

_backend_upload_pairs!(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc} =
    (cff.gpu isa GPUState && _upload_pairs!(cff); nothing)

function _upload_positions!(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc}
    g = cff.gpu::GPUState
    @inbounds for k in 1:cff.natoms
        rk = cff.r[k]; g.hx[k] = rk[1]; g.hy[k] = rk[2]; g.hz[k] = rk[3]
    end
    copyto!(g.rx, g.hx); copyto!(g.ry, g.hy); copyto!(g.rz, g.hz)
    nothing
end

# ---------------------------------------------------------------------------
#  Kernels (one source; per-pair math shared with the CPU loops)
# ---------------------------------------------------------------------------
@kernel function _k_force_lj!(Fx, Fy, Fz, @Const(rx), @Const(ry), @Const(rz),
        @Const(I), @Const(J), @Const(A), @Const(B), @Const(S), sw, N)
    n = @index(Global)
    if n <= N
        @inbounds begin
            i = I[n]; j = J[n]
            dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j]
            factor = _f_lj_factor(dx*dx + dy*dy + dz*dz, A[n], B[n], S[n], sw)
            fx = factor*dx; fy = factor*dy; fz = factor*dz
            Atomix.@atomic Fx[i] += fx; Atomix.@atomic Fy[i] += fy; Atomix.@atomic Fz[i] += fz
            Atomix.@atomic Fx[j] -= fx; Atomix.@atomic Fy[j] -= fy; Atomix.@atomic Fz[j] -= fz
        end
    end
end

@kernel function _k_force_hb!(Fx, Fy, Fz, @Const(rx), @Const(ry), @Const(rz),
        @Const(I), @Const(J), @Const(A), @Const(B), @Const(S), sw, N)
    n = @index(Global)
    if n <= N
        @inbounds begin
            i = I[n]; j = J[n]
            dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j]
            factor = _f_hb_factor(dx*dx + dy*dy + dz*dz, A[n], B[n], S[n], sw)
            fx = factor*dx; fy = factor*dy; fz = factor*dz
            Atomix.@atomic Fx[i] += fx; Atomix.@atomic Fy[i] += fy; Atomix.@atomic Fz[i] += fz
            Atomix.@atomic Fx[j] -= fx; Atomix.@atomic Fy[j] -= fy; Atomix.@atomic Fz[j] -= fz
        end
    end
end

@kernel function _k_force_es!(Fx, Fy, Fz, @Const(rx), @Const(ry), @Const(rz),
        @Const(I), @Const(J), @Const(Q), @Const(S), sw, ddd, pref, N)
    n = @index(Global)
    if n <= N
        @inbounds begin
            i = I[n]; j = J[n]
            dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j]
            factor = _f_es_factor(dx*dx + dy*dy + dz*dz, Q[n], S[n], sw, ddd, pref)
            fx = factor*dx; fy = factor*dy; fz = factor*dz
            Atomix.@atomic Fx[i] += fx; Atomix.@atomic Fy[i] += fy; Atomix.@atomic Fz[i] += fz
            Atomix.@atomic Fx[j] -= fx; Atomix.@atomic Fy[j] -= fy; Atomix.@atomic Fz[j] -= fz
        end
    end
end

@kernel function _k_energy_lj!(E, @Const(rx), @Const(ry), @Const(rz),
        @Const(I), @Const(J), @Const(A), @Const(B), @Const(S), sw, N)
    n = @index(Global)
    if n <= N
        @inbounds begin
            i = I[n]; j = J[n]
            dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j]
            E[n] = _e_lj_pair(dx*dx + dy*dy + dz*dz, A[n], B[n], S[n], sw)
        end
    end
end

@kernel function _k_energy_hb!(E, @Const(rx), @Const(ry), @Const(rz),
        @Const(I), @Const(J), @Const(A), @Const(B), @Const(S), sw, N)
    n = @index(Global)
    if n <= N
        @inbounds begin
            i = I[n]; j = J[n]
            dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j]
            E[n] = _e_hb_pair(dx*dx + dy*dy + dz*dz, A[n], B[n], S[n], sw)
        end
    end
end

@kernel function _k_energy_es!(E, @Const(rx), @Const(ry), @Const(rz),
        @Const(I), @Const(J), @Const(Q), @Const(S), sw, ddd, pref, N)
    n = @index(Global)
    if n <= N
        @inbounds begin
            i = I[n]; j = J[n]
            dx = rx[i] - rx[j]; dy = ry[i] - ry[j]; dz = rz[i] - rz[j]
            E[n] = _e_es_pair(dx*dx + dy*dy + dz*dz, Q[n], S[n], sw, ddd, pref)
        end
    end
end

# ---------------------------------------------------------------------------
#  GPU evaluation entry points
# ---------------------------------------------------------------------------
function compute_energy!(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc}
    dev = _dev(cff); wg = _wg(cff); g = cff.gpu::GPUState
    _upload_positions!(cff)
    nlj = length(cff.lj.i); nhb = length(cff.hb.i); nes = length(cff.es.i)
    nlj > 0 && _k_energy_lj!(dev, wg)(g.lj_e, g.rx, g.ry, g.rz, g.lj_i, g.lj_j, g.lj_A, g.lj_B, g.lj_s, cff.vdw_sw, nlj; ndrange=nlj)
    nhb > 0 && _k_energy_hb!(dev, wg)(g.hb_e, g.rx, g.ry, g.rz, g.hb_i, g.hb_j, g.hb_A, g.hb_B, g.hb_s, cff.vdw_sw, nhb; ndrange=nhb)
    nes > 0 && _k_energy_es!(dev, wg)(g.es_e, g.rx, g.ry, g.rz, g.es_i, g.es_j, g.es_q, g.es_s, cff.es_sw, cff.distance_dependent_dielectric, cff.es_prefactor, nes; ndrange=nes)
    _KA.synchronize(dev)
    vdw   = nlj > 0 ? Acc(sum(g.lj_e)) : zero(Acc)
    hbond = nhb > 0 ? Acc(sum(g.hb_e)) : zero(Acc)
    es    = nes > 0 ? Acc(sum(g.es_e)) : zero(Acc)
    _finish_energy!(cff, vdw, hbond, es)
end

function compute_forces!(cff::CompiledForceField{T, Acc, <:KAGPUBackend}) where {T, Acc}
    dev = _dev(cff); wg = _wg(cff); g = cff.gpu::GPUState
    _upload_positions!(cff)
    fill!(g.Fx, zero(Acc)); fill!(g.Fy, zero(Acc)); fill!(g.Fz, zero(Acc))
    nlj = length(cff.lj.i); nhb = length(cff.hb.i); nes = length(cff.es.i)
    nlj > 0 && _k_force_lj!(dev, wg)(g.Fx, g.Fy, g.Fz, g.rx, g.ry, g.rz, g.lj_i, g.lj_j, g.lj_A, g.lj_B, g.lj_s, cff.vdw_sw, nlj; ndrange=nlj)
    nhb > 0 && _k_force_hb!(dev, wg)(g.Fx, g.Fy, g.Fz, g.rx, g.ry, g.rz, g.hb_i, g.hb_j, g.hb_A, g.hb_B, g.hb_s, cff.vdw_sw, nhb; ndrange=nhb)
    nes > 0 && _k_force_es!(dev, wg)(g.Fx, g.Fy, g.Fz, g.rx, g.ry, g.rz, g.es_i, g.es_j, g.es_q, g.es_s, cff.es_sw, cff.distance_dependent_dielectric, cff.es_prefactor_force, nes; ndrange=nes)
    _KA.synchronize(dev)
    copyto!(g.hFx, g.Fx); copyto!(g.hFy, g.Fy); copyto!(g.hFz, g.Fz)
    F = cff.F
    @inbounds for k in 1:cff.natoms
        F[k] = Vector3{Acc}(g.hFx[k], g.hFy[k], g.hFz[k])
    end
    _force_bonded!(F, cff)     # cheap, irregular -> host
    _apply_constraints!(cff)
    nothing
end

end # module
