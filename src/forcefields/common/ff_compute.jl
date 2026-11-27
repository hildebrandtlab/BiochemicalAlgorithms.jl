@inline function _compute_stretch_energy(k::T, r0::T, distance::T)::T where T
    k * (distance - r0)^2
end

@inline function _compute_stretch_force(k::T, r0::T, distance::T, direction::Vector3{T})::Vector3{T} where T
    if distance == zero(T)
        return zero(direction)
    end
    2 * k * (distance - r0) / distance * direction
end

@inline function _compute_bend_energy(k::T, θ0::T, v1::Vector3{T}, v2::Vector3{T})::T where T
    sq_length = squared_norm(v1) * squared_norm(v2)

    if iszero(sq_length)
        return zero(T)
    end

    cos_θ = dot(v1, v2) / sqrt(sq_length)
    θ = if cos_θ > one(T)
            zero(T)
        elseif cos_θ < -one(T)
            π
        else
            acos(cos_θ)
        end

    k * (θ - θ0)^2
end

@inline function _compute_bend_force(k::T, θ0::T, v1::Vector3{T}, v2::Vector3{T})::Tuple{Vector3{T}, Vector3{T}, Vector3{T}} where T
    v1_length = norm(v1)
    v2_length = norm(v2)

    if v1_length == zero(T) || v2_length == zero(T)
        return zero(v1), zero(v1), zero(v1)
    end

    v1_norm = v1 / v1_length
    v2_norm = v2 / v2_length

    cos_θ = dot(v1_norm, v2_norm)
    θ = if cos_θ > one(T)
            zero(T)
        elseif cos_θ < -one(T)
            π
        else
            acos(cos_θ)
        end

    factor = T(2) * k * (θ - θ0)

    crossv1v2 = normalize(cross(v1_norm, v2_norm))

    if isnan(crossv1v2[1])
        return zero(v1), zero(v1), zero(v1)
    end

    n1 = (cross(v1_norm, crossv1v2)) .* factor ./ v1_length
    n2 = (cross(v2_norm, crossv1v2)) .* factor ./ v2_length

    f1 = -n1
    f2 = n1 - n2
    f3 = n2

    f1, f2, f3
end

@inline function _torsion_geometry(a21::Vector3{T}, a23::Vector3{T}, a34::Vector3{T}) where T
    cross2321 = cross(a23, a21)
    cross2334 = cross(a23, a34)
    len1 = norm(cross2321); len2 = norm(cross2334)
    ok = !(len1 == zero(T) || len2 == zero(T))
    ncross2321 = ok ? cross2321 / len1 : cross2321
    ncross2334 = ok ? cross2334 / len2 : cross2334
    cos_ϕ = ok ? clamp(dot(ncross2321, ncross2334), -one(T), one(T)) : zero(T)
    (cross2321, cross2334, len1, len2, ok, cos_ϕ)
end

@inline function _compute_torsion_energy(V::AbstractVector{T}, ϕ₀::AbstractVector{T}, f::AbstractVector{<:Integer}, div::AbstractVector{<:Integer}, a21::Vector3{T}, a23::Vector3{T}, a34::Vector3{T})::T where T
    _, _, _, _, ok, cos_ϕ = _torsion_geometry(a21, a23, a34)
    ok || return zero(T)

    acϕ = acos(cos_ϕ)
    energy = zero(T)
    @inbounds for n in eachindex(V)
        energy += V[n] / T(div[n]) * (one(T) + cos(T(f[n]) * acϕ - ϕ₀[n]))
    end
    energy
end

@inline function _compute_torsion_force(V::AbstractVector{T}, ϕ₀::AbstractVector{T}, f::AbstractVector{<:Integer}, div::AbstractVector{<:Integer}, a21::Vector3{T}, a23::Vector3{T}, a34::Vector3{T})::Tuple{Vector3{T}, Vector3{T}, Vector3{T}, Vector3{T}} where T
    cross2321, cross2334, len1, len2, ok, cos_ϕ = _torsion_geometry(a21, a23, a34)
    ok || return zero(a23), zero(a23), zero(a23), zero(a23)

    acϕ = acos(cos_ϕ)
    ∂E∂ϕ = zero(T)
    @inbounds for n in eachindex(V)
        ∂E∂ϕ -= V[n] / T(div[n]) * T(f[n]) * sin(T(f[n]) * acϕ - ϕ₀[n])
    end

    direction = dot(cross(cross2321, cross2334), a23)
    direction > zero(T) && (∂E∂ϕ *= -one(T))

    a13 = a23 .- a21  # a3 - a1 = (a3 - a2) - (a2 - a1)
    a24 = a34 .+ a23  # a4 - a2 = (a4 - a3) + (a3 - a2)
    na23 = norm(a23)

    dEdt =  (∂E∂ϕ / (len1^2 * na23)) * cross(cross2321, a23)
    dEdu = -(∂E∂ϕ / (len2^2 * na23)) * cross(cross2334, a23)

    f1 = cross(dEdt, a23)
    f2 = cross(a13, dEdt) + cross(dEdu, a34)
    f3 = cross(a21, dEdt) + cross(a24, dEdu)
    f4 = cross(dEdu, a23)

    f1, f2, f3, f4
end

# ---------------------------------------------------------------------------
#  Per-pair nonbonded math (scalar; shared by reference / serial / threaded / GPU).
#
#  Energy helpers return the pair's energy contribution given the squared
#  distance `sq`. Force helpers return a scalar `factor` such that the force on
#  atom i is `factor * (r_i - r_j)` (and `-factor * (r_i - r_j)` on j), so a
#  single scalar unifies the CPU (SVector) and GPU (scalar SoA) loops.
# ---------------------------------------------------------------------------
@inline function _e_lj_pair(sq::Acc, A, B, scaling, sw::CubicSwitchingFunction{Acc}) where Acc
    inv6 = sq^-3
    inv6 * (inv6 * A - B) * scaling * switching_function(sw, sq)
end

@inline function _e_hb_pair(sq::Acc, A, B, scaling, sw::CubicSwitchingFunction{Acc}) where Acc
    (sq^-6 * A - sq^-5 * B) * scaling * switching_function(sw, sq)
end

@inline function _e_es_pair(sq::Acc, q1q2, scaling, sw::CubicSwitchingFunction{Acc},
                            ddd::Bool, pref::Acc) where Acc
    inner = ddd ? q1q2 / 4 / sq : q1q2 / sqrt(sq)
    inner * scaling * switching_function(sw, sq) * pref
end

@inline function _f_lj_factor(sq::Acc, A, B, scaling, sw::CubicSwitchingFunction{Acc}) where Acc
    (sq > zero(Acc) && sq <= sw.sq_cutoff) || return zero(Acc)
    inv6 = sq^-3
    factor = (one(Acc) / sq) * inv6 * scaling * (12 * A * inv6 - 6 * B)
    if sq > sw.sq_cuton
        sval, sder = switching_derivative(sw, sq)
        factor *= sval
        factor += sder * (-scaling * inv6 * (inv6 * A - B))
    end
    factor
end

@inline function _f_hb_factor(sq::Acc, A, B, scaling, sw::CubicSwitchingFunction{Acc}) where Acc
    (sq > zero(Acc) && sq <= sw.sq_cutoff) || return zero(Acc)
    inv2 = one(Acc) / sq
    factor = inv2 * sq^-5 * scaling * (12 * A * inv2 - 10 * B)
    if sq > sw.sq_cuton
        sval, sder = switching_derivative(sw, sq)
        factor *= sval
        factor += sder * (-scaling * sq^-5 * (A * inv2 - B))
    end
    factor
end

@inline function _f_es_factor(sq::Acc, q1q2, scaling, sw::CubicSwitchingFunction{Acc},
                              ddd::Bool, pref::Acc) where Acc
    (sq > zero(Acc) && sq <= sw.sq_cutoff) || return zero(Acc)
    inv2 = one(Acc) / sq
    inv = sqrt(inv2)
    factor = q1q2 * inv2 * scaling * pref
    factor *= ddd ? Acc(0.5) * inv2 : inv
    if sq > sw.sq_cuton
        sval, sder = switching_derivative(sw, sq)
        factor *= sval
        ddf = ddd ? Acc(0.25) * inv : one(Acc)
        factor += sder * (-pref * scaling * ddf * inv * q1q2)
    end
    factor
end

# Minimum-image displacement for an orthorhombic box of side lengths `L`
# (valid for cutoff < L/2). `L[1] <= 0` marks a non-periodic system.
@inline function _disp(a::Vector3{Acc}, b::Vector3{Acc}, L::Vector3{Acc}) where Acc
    d = a - b
    L[1] > zero(Acc) ? d - L .* round.(d ./ L) : d
end

@inline function _compute_lj_energy(A::T, B::T, distance::T, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::T where T
    _e_lj_pair(distance^2, A, B, scaling_factor, switch_func)
end

@inline function _compute_hbond_energy(A::T, B::T, distance::T, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::T where T
    _e_hb_pair(distance^2, A, B, scaling_factor, switch_func)
end

@inline function _compute_es_energy(q1q2::T, distance::T, scaling_factor::T, distance_dependent_dielectric::Bool, switch_func::CubicSwitchingFunction{T})::T where T
    _e_es_pair(distance^2, q1q2, scaling_factor, switch_func, distance_dependent_dielectric, T(ES_Prefactor))
end

@inline function _compute_lj_force(A::T, B::T, sq_distance::T, direction::Vector3{T}, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::Vector3{T} where T
    _f_lj_factor(sq_distance, A, B, scaling_factor, switch_func) * direction
end

@inline function _compute_hbond_force(A::T, B::T, sq_distance::T, direction::Vector3{T}, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::Vector3{T} where T
    _f_hb_factor(sq_distance, A, B, scaling_factor, switch_func) * direction
end

@inline function _compute_es_force(q1q2::T, sq_distance::T, direction::Vector3{T}, scaling_factor::T, distance_dependent_dielectric::Bool, switch_func::CubicSwitchingFunction{T})::Vector3{T} where T
    _f_es_factor(sq_distance, q1q2, scaling_factor, switch_func, distance_dependent_dielectric, T(ES_Prefactor_force)) * direction
end
