export
    _compute_stretch_energy,
    _compute_stretch_force,
    _compute_bend_energy,
    _compute_bend_force,
    _torsion_geometry,
    _compute_torsion_energy,
    _compute_torsion_force,
    _compute_lj_energy,
    _compute_hbond_energy,
    _compute_es_energy,
    _compute_lj_force,
    _compute_hbond_force,
    _compute_es_force

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

@inline function _compute_torsion_energy(V::Union{T, AbstractVector{T}}, ϕ₀::Union{T, AbstractVector{T}}, f::Union{Integer, AbstractVector{<:Integer}}, div::Union{Integer, AbstractVector{<:Integer}}, a21::Vector3{T}, a23::Vector3{T}, a34::Vector3{T})::T where T
    _, _, _, _, ok, cos_ϕ = _torsion_geometry(a21, a23, a34)
    ok || return zero(T)

    acϕ = acos(cos_ϕ)
    if V isa AbstractVector
        terms = V ./ T.(div) .* (Ref(1) .+ cos.(T.(f) .* Ref(acϕ) .- ϕ₀))
        return sum(terms)
    else
        return V / T(div) * (1 + cos(T(f) * acϕ - ϕ₀))
    end
end

@inline function _compute_torsion_force(V::Union{T, AbstractVector{T}}, ϕ₀::Union{T, AbstractVector{T}}, f::Union{Integer, AbstractVector{<:Integer}}, div::Union{Integer, AbstractVector{<:Integer}}, a21::Vector3{T}, a23::Vector3{T}, a34::Vector3{T})::Tuple{Vector3{T}, Vector3{T}, Vector3{T}, Vector3{T}} where T
    cross2321, cross2334, len1, len2, ok, cos_ϕ = _torsion_geometry(a21, a23, a34)
    ok || return zero(a23), zero(a23), zero(a23), zero(a23)

    if V isa AbstractVector
        terms = -V./T.(div) .* T.(f) .* (sin.(T.(f) .* Ref(acos(cos_ϕ)) .- ϕ₀))
        ∂E∂ϕ = sum(terms)
    else
        ∂E∂ϕ = -V / T(div) * T(f) * sin(T(f) * acos(cos_ϕ) - ϕ₀)
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

@inline function _compute_lj_energy(A::T, B::T, distance::T, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::T where T
    inv_dist_6::T = distance^-6
    switch_value = switching_function(switch_func, distance^2)
    scaling_factor * inv_dist_6 * (inv_dist_6 * A - B) * switch_value
end

@inline function _compute_hbond_energy(A::T, B::T, distance::T, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::T where T
    switch_value = switching_function(switch_func, distance^2)
    scaling_factor * (distance^-12 * A - distance^-10 * B) * switch_value
end

@inline function _compute_es_energy(q1q2::T, distance::T, scaling_factor::T, distance_dependent_dielectric::Bool, switch_func::CubicSwitchingFunction{T})::T where T
    energy::T = distance_dependent_dielectric ? q1q2 / 4 / distance^2 : q1q2 / distance
    energy * scaling_factor * switching_function(switch_func, distance^2) * T(ES_Prefactor)
end

@inline function _compute_lj_force(A::T, B::T, sq_distance::T, direction::Vector3{T}, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::Vector3{T} where T
    if (sq_distance > zero(T) && sq_distance <= switch_func.sq_cutoff)
        inv_distance_2 = one(T) / sq_distance
        inv_distance_6 = sq_distance^-3

        factor = inv_distance_2 * inv_distance_6 * scaling_factor * (12 * A * inv_distance_6 - 6 * B)

        if (sq_distance > switch_func.sq_cuton)
            switch_value, switch_derivative = switching_derivative(switch_func, sq_distance)
            factor *= switch_value

            energy = -scaling_factor * inv_distance_6 * (inv_distance_6 * A - B)
            factor += switch_derivative * energy
        end

        return factor * direction
    end
    return zero(direction)
end

@inline function _compute_hbond_force(A::T, B::T, sq_distance::T, direction::Vector3{T}, scaling_factor::T, switch_func::CubicSwitchingFunction{T})::Vector3{T} where T
    if (sq_distance > zero(T) && sq_distance <= switch_func.sq_cutoff)
        inv_distance_2 = one(T) / sq_distance
        inv_distance_10 = sq_distance^-5
        inv_distance_12 = sq_distance^-6

        factor = inv_distance_2 * inv_distance_12 * (12 * A * inv_distance_2 - 10 * B) * scaling_factor

        if (sq_distance > switch_func.sq_cuton)
            switch_value, switch_derivative = switching_derivative(switch_func, sq_distance)
            factor *= switch_value

            energy = -scaling_factor * inv_distance_10 * (A * inv_distance_2 - B)
            factor += switch_derivative * energy
        end

        return factor * direction
    end
    return zero(direction)
end

@inline function _compute_es_force(q1q2::T, sq_distance::T, direction::Vector3{T}, scaling_factor::T, distance_dependent_dielectric::Bool, switch_func::CubicSwitchingFunction{T})::Vector3{T} where T
    if (sq_distance > zero(T) && sq_distance <= switch_func.sq_cutoff)
        inv_distance_2 = one(T) / sq_distance
        inv_distance = sqrt(inv_distance_2)

        factor = q1q2 * inv_distance_2 * scaling_factor * T(ES_Prefactor_force)

        if distance_dependent_dielectric
            factor *= T(0.5) * inv_distance_2
        else
            factor *= inv_distance
        end

        if (sq_distance > switch_func.sq_cuton)
            switch_value, switch_derivative = switching_derivative(switch_func, sq_distance)
            factor *= switch_value

            dist_depend_factor = distance_dependent_dielectric ? T(0.25) * inv_distance : one(T)
            energy = -T(ES_Prefactor_force) * scaling_factor * dist_depend_factor * inv_distance * q1q2
            factor += switch_derivative * energy
        end

        return factor * direction
    end
    return zero(direction)
end
