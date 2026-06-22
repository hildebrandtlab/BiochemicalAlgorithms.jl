# ============================================================================
#  Energy / force kernels for the CompiledForceField.
#
#  Every formula here is a verbatim transcription of the corresponding
#  `compute_energy` / `compute_forces!` in the reference components, rewritten
#  to read positions from `cff.r[row]` and accumulate forces into `cff.F[row]`
#  (no `Atom{T}` row-views, no Dict lookups). With acc == storage type the
#  results match the reference path bit-for-bit.
#
#  Stage 1 provides the SerialBackend; later stages add threaded / GPU loops
#  that call the same per-interaction math.
# ============================================================================

# ---------------------------------------------------------------------------
#  Per-component ENERGY (serial accumulation in Acc)
# ---------------------------------------------------------------------------
function _energy_stretch(cff::CompiledForceField{T, Acc}) where {T, Acc}
    s = cff.stretch; r = cff.r
    e = zero(Acc)
    @inbounds for n in eachindex(s.k)
        d = norm(r[s.i[n]] - r[s.j[n]])
        e += s.k[n] * (d - s.r0[n])^2
    end
    e
end

function _energy_bend(cff::CompiledForceField{T, Acc}) where {T, Acc}
    b = cff.bend; r = cff.r
    e = zero(Acc)
    @inbounds for n in eachindex(b.k)
        v1 = r[b.i[n]] - r[b.j[n]]
        v2 = r[b.l[n]] - r[b.j[n]]
        sq_length = squared_norm(v1) * squared_norm(v2)
        iszero(sq_length) && continue
        cos_θ = dot(v1, v2) / sqrt(sq_length)
        θ = cos_θ > one(Acc) ? zero(Acc) : cos_θ < -one(Acc) ? Acc(π) : acos(cos_θ)
        e += b.k[n] * (θ - b.θ0[n])^2
    end
    e
end

function _energy_torsion(cff::CompiledForceField{T, Acc}, t::TorsionArrays{Acc}) where {T, Acc}
    r = cff.r
    e = zero(Acc)
    @inbounds for n in eachindex(t.i)
        a23 = r[t.k[n]] - r[t.j[n]]
        cross2321 = normalize(cross(a23, r[t.i[n]] - r[t.j[n]]))
        cross2334 = normalize(cross(a23, r[t.l[n]] - r[t.k[n]]))
        (isnan(cross2321[1]) || isnan(cross2334[1])) && continue
        cos_ϕ = clamp(dot(cross2321, cross2334), -one(Acc), one(Acc))
        acϕ = acos(cos_ϕ)
        lo = t.term_offset[n]; hi = t.term_offset[n+1] - one(Int32)
        for m in lo:hi
            e += t.V[m] / t.div[m] * (one(Acc) + cos(t.f[m] * acϕ - t.ϕ0[m]))
        end
    end
    e
end

function _energy_lj(cff::CompiledForceField{T, Acc}) where {T, Acc}
    p = cff.lj; r = cff.r; sw = cff.vdw_sw
    e = zero(Acc)
    @inbounds for n in eachindex(p.i)
        sq = squared_norm(r[p.i[n]] - r[p.j[n]])
        inv6 = sq^-3
        e += inv6 * (inv6 * p.A[n] - p.B[n]) * p.scaling[n] * switching_function(sw, sq)
    end
    e
end

function _energy_hb(cff::CompiledForceField{T, Acc}) where {T, Acc}
    p = cff.hb; r = cff.r; sw = cff.vdw_sw
    e = zero(Acc)
    @inbounds for n in eachindex(p.i)
        sq = squared_norm(r[p.i[n]] - r[p.j[n]])
        # NOTE: precedence transcribed verbatim from the reference (scaling and
        # switching multiply only the second term).
        e += sq^-6 * p.A[n] - sq^-5 * p.B[n] * p.scaling[n] * switching_function(sw, sq)
    end
    e
end

function _energy_es(cff::CompiledForceField{T, Acc}) where {T, Acc}
    p = cff.es; r = cff.r; sw = cff.es_sw
    ddd = cff.distance_dependent_dielectric
    pref = cff.es_prefactor
    e = zero(Acc)
    @inbounds for n in eachindex(p.i)
        sq = squared_norm(r[p.i[n]] - r[p.j[n]])
        inner = ddd ? p.q1q2[n] / 4 / sq : p.q1q2[n] / sqrt(sq)
        e += inner * p.scaling[n] * switching_function(sw, sq) * pref
    end
    e
end

# ---------------------------------------------------------------------------
#  Per-component FORCES (serial accumulation into cff.F)
# ---------------------------------------------------------------------------
function _force_stretch!(cff::CompiledForceField{T, Acc}) where {T, Acc}
    s = cff.stretch; r = cff.r; F = cff.F
    @inbounds for n in eachindex(s.k)
        i = s.i[n]; j = s.j[n]
        direction = r[i] - r[j]
        d = norm(direction)
        d == zero(Acc) && continue
        direction *= 2 * s.k[n] * (d - s.r0[n]) / d
        F[i] -= direction
        F[j] += direction
    end
    nothing
end

function _force_bend!(cff::CompiledForceField{T, Acc}) where {T, Acc}
    b = cff.bend; r = cff.r; F = cff.F
    @inbounds for n in eachindex(b.k)
        i = b.i[n]; j = b.j[n]; l = b.l[n]
        v1 = r[i] - r[j]
        v2 = r[l] - r[j]
        v1_length = norm(v1); v2_length = norm(v2)
        (v1_length == zero(Acc) || v2_length == zero(Acc)) && continue
        v1 = v1 / v1_length
        v2 = v2 / v2_length
        cos_θ = dot(v1, v2)
        θ = cos_θ > one(Acc) ? zero(Acc) : cos_θ < -one(Acc) ? Acc(π) : acos(cos_θ)
        factor = 2 * b.k[n] * (θ - b.θ0[n])
        crossv1v2 = normalize(cross(v1, v2))
        isnan(crossv1v2[1]) && continue
        n1 = cross(v1, crossv1v2) .* factor ./ v1_length
        n2 = cross(v2, crossv1v2) .* factor ./ v2_length
        F[i] -= n1
        F[j] += n1
        F[j] -= n2
        F[l] += n2
    end
    nothing
end

function _force_torsion!(cff::CompiledForceField{T, Acc}, t::TorsionArrays{Acc}) where {T, Acc}
    r = cff.r; F = cff.F
    @inbounds for n in eachindex(t.i)
        i = t.i[n]; j = t.j[n]; k = t.k[n]; l = t.l[n]
        a21 = r[i] - r[j]
        a23 = r[k] - r[j]
        a34 = r[l] - r[k]
        cross2321 = cross(a23, a21)
        cross2334 = cross(a23, a34)
        len1 = norm(cross2321); len2 = norm(cross2334)
        (len1 == zero(Acc) || len2 == zero(Acc)) && continue
        cos_ϕ = clamp(dot(cross2321, cross2334) / (len1 * len2), -one(Acc), one(Acc))
        acϕ = acos(cos_ϕ)
        lo = t.term_offset[n]; hi = t.term_offset[n+1] - one(Int32)
        ∂E∂ϕ = zero(Acc)
        for m in lo:hi
            ∂E∂ϕ += -t.V[m] / t.div[m] * t.f[m] * sin(t.f[m] * acϕ - t.ϕ0[m])
        end
        direction = dot(cross(cross2321, cross2334), a23)
        direction > zero(Acc) && (∂E∂ϕ = -∂E∂ϕ)
        a13 = r[k] - r[i]
        a24 = r[l] - r[j]
        na23 = norm(a23)
        dEdt =  (∂E∂ϕ / (len1^2 * na23)) * cross(cross2321, a23)
        dEdu = -(∂E∂ϕ / (len2^2 * na23)) * cross(cross2334, a23)
        F[i] += cross(dEdt, a23)
        F[j] += cross(a13, dEdt) + cross(dEdu, a34)
        F[k] += cross(a21, dEdt) + cross(a24, dEdu)
        F[l] += cross(dEdu, a23)
    end
    nothing
end

function _force_lj!(cff::CompiledForceField{T, Acc}) where {T, Acc}
    p = cff.lj; r = cff.r; F = cff.F; sw = cff.vdw_sw
    @inbounds for n in eachindex(p.i)
        i = p.i[n]; j = p.j[n]
        direction = r[i] - r[j]
        sq = squared_norm(direction)
        (sq > zero(Acc) && sq <= sw.sq_cutoff) || continue
        factor = one(Acc) / sq
        inv6 = sq^-3
        factor *= inv6 * p.scaling[n] * (12 * p.A[n] * inv6 - 6 * p.B[n])
        if sq > sw.sq_cuton
            sval, sder = switching_derivative(sw, sq)
            factor *= sval
            energy = -p.scaling[n] * inv6 * (inv6 * p.A[n] - p.B[n])
            factor += sder * energy
        end
        force = factor * direction
        F[i] += force
        F[j] -= force
    end
    nothing
end

function _force_hb!(cff::CompiledForceField{T, Acc}) where {T, Acc}
    p = cff.hb; r = cff.r; F = cff.F; sw = cff.vdw_sw
    @inbounds for n in eachindex(p.i)
        i = p.i[n]; j = p.j[n]
        direction = r[i] - r[j]
        sq = squared_norm(direction)
        (sq > zero(Acc) && sq <= sw.sq_cutoff) || continue
        inv2 = one(Acc) / sq
        inv10 = sq^-5
        inv12 = sq^-6
        factor = inv2
        factor *= inv12 * (12 * p.A[n] * inv2 - 10 * p.B[n])
        if sq > sw.sq_cuton
            sval, sder = switching_derivative(sw, sq)
            factor *= sval
            energy = -p.scaling[n] * inv10 * (p.A[n] * inv2 - p.B[n])
            factor += sder * energy
        end
        force = factor * direction
        F[i] += force
        F[j] -= force
    end
    nothing
end

function _force_es!(cff::CompiledForceField{T, Acc}) where {T, Acc}
    p = cff.es; r = cff.r; F = cff.F; sw = cff.es_sw
    ddd = cff.distance_dependent_dielectric
    pref = cff.es_prefactor_force
    @inbounds for n in eachindex(p.i)
        i = p.i[n]; j = p.j[n]
        direction = r[i] - r[j]
        sq = squared_norm(direction)
        (sq > zero(Acc) && sq <= sw.sq_cutoff) || continue
        inv2 = one(Acc) / sq
        inv = sqrt(inv2)
        factor = p.q1q2[n] * inv2 * p.scaling[n] * pref
        if ddd
            factor *= Acc(0.5) * inv2
        else
            factor *= inv
        end
        if sq > sw.sq_cuton
            sval, sder = switching_derivative(sw, sq)
            factor *= sval
            ddf = ddd ? Acc(0.25) * inv : one(Acc)
            energy = -pref * p.scaling[n] * ddf * inv * p.q1q2[n]
            factor += sder * energy
        end
        force = factor * direction
        F[i] += force
        F[j] -= force
    end
    nothing
end

# ---------------------------------------------------------------------------
#  Public evaluation entry points (dispatch on backend)
# ---------------------------------------------------------------------------
"""
    compute_energy!(cff::CompiledForceField) -> Acc

Evaluate the total energy at the evaluator's current coordinates, filling
`cff.energy` with the per-category breakdown (same keys as the reference
`ForceField` path).
"""
function compute_energy!(cff::CompiledForceField{T, Acc, SerialBackend}) where {T, Acc}
    stretch   = _energy_stretch(cff)
    bend      = _energy_bend(cff)
    proper    = _energy_torsion(cff, cff.proper)
    improper  = _energy_torsion(cff, cff.improper)
    vdw       = _energy_lj(cff)
    hbond     = _energy_hb(cff)
    es        = _energy_es(cff)

    cff.energy["Bond Stretches"]   = stretch
    cff.energy["Angle Bends"]      = bend
    cff.energy["Proper Torsion"]   = proper
    cff.energy["Improper Torsion"] = improper
    cff.energy["Van der Waals"]    = vdw
    cff.energy["Hydrogen Bonds"]   = hbond
    cff.energy["Electrostatic"]    = es

    stretch + bend + proper + improper + vdw + hbond + es
end

"""
    compute_forces!(cff::CompiledForceField)

Zero and recompute forces at the current coordinates, accumulating into
`cff.F`. Forces on constrained atoms are zeroed afterwards.
"""
function compute_forces!(cff::CompiledForceField{T, Acc, SerialBackend}) where {T, Acc}
    fill!(cff.F, zero(Vector3{Acc}))
    _force_stretch!(cff)
    _force_bend!(cff)
    _force_torsion!(cff, cff.proper)
    _force_torsion!(cff, cff.improper)
    _force_lj!(cff)
    _force_hb!(cff)
    _force_es!(cff)
    _apply_constraints!(cff)
    nothing
end

@inline function _apply_constraints!(cff::CompiledForceField{T, Acc}) where {T, Acc}
    any(cff.constrained_mask) || return nothing
    @inbounds for k in 1:cff.natoms
        cff.constrained_mask[k] && (cff.F[k] = zero(Vector3{Acc}))
    end
    nothing
end
