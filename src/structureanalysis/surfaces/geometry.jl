"""
3D analytical geometry primitives used by the Connolly-style molecular surface
algorithms. All routines return `nothing` for the degenerate cases (e.g. no
intersection) and propagate the type parameter `T`.

These mirror BALL's `BALL/MATHS/analyticalGeometry.h` for the entry points the
surface algorithms actually call.
"""

# ---------------------------------------------------------------------------
# Sphere–sphere intersection: returns the circle along which the two spheres
# touch, or `nothing` if they don't intersect.
#
# The plane of intersection is perpendicular to (c2 - c1). Let d = |c2 - c1|.
# The signed distance from c1 to the plane is
#     a = (d² + r1² − r2²) / (2 d)
# and the radius of the circle is sqrt(r1² − a²).
# ---------------------------------------------------------------------------
@inline function intersect_spheres(s1::Sphere{T}, s2::Sphere{T}) where {T<:Real}
    d_vec = s2.center - s1.center
    d2 = dot(d_vec, d_vec)
    iszero(d2) && return nothing
    d = sqrt(d2)
    a = (d2 + s1.r^2 - s2.r^2) / (2 * d)
    h2 = s1.r^2 - a^2
    h2 < 0 && return nothing
    n = d_vec / d
    p = s1.center + a * n
    Circle3{T}(Vector3{T}(p[1], p[2], p[3]), Vector3{T}(n[1], n[2], n[3]), sqrt(h2))
end

# ---------------------------------------------------------------------------
# Sphere–sphere–sphere intersection: returns the two candidate points where all
# three spheres meet (or `nothing` if there is no real intersection).
#
# Work in a local frame with c1 at the origin, c2 on the +x axis, c3 in the xy
# plane. Solve the 3×3 system in closed form, then transform back. The two
# solutions differ in the sign of the z component (above / below the plane).
# ---------------------------------------------------------------------------
function intersect_three_spheres(s1::Sphere{T}, s2::Sphere{T}, s3::Sphere{T}) where {T<:Real}
    e1 = s2.center - s1.center
    d = norm(e1)
    iszero(d) && return nothing
    ex = e1 / d

    e2 = s3.center - s1.center
    i = dot(ex, e2)
    e2_perp = e2 - i * ex
    j = norm(e2_perp)
    iszero(j) && return nothing
    ey = e2_perp / j

    ez = cross(ex, ey)

    r1sq = s1.r^2
    x = (r1sq - s2.r^2 + d*d) / (2 * d)
    y = (r1sq - s3.r^2 + i*i + j*j - 2 * i * x) / (2 * j)
    z2 = r1sq - x*x - y*y
    z2 < 0 && return nothing
    z = sqrt(z2)

    base = s1.center + x * ex + y * ey
    off  = z * ez
    p1 = base + off
    p2 = base - off
    (Vector3{T}(p1[1], p1[2], p1[3]), Vector3{T}(p2[1], p2[2], p2[3]))
end

# ---------------------------------------------------------------------------
# Sphere–line intersection: parametric line `o + t·d` (d need not be unit).
# Returns the two intersection points (ordered by increasing `t`) or `nothing`.
# ---------------------------------------------------------------------------
function intersect_sphere_line(s::Sphere{T}, o::AbstractVector{T}, d::AbstractVector{T}) where {T<:Real}
    f = o - s.center
    a = dot(d, d)
    iszero(a) && return nothing
    b = 2 * dot(f, d)
    c = dot(f, f) - s.r^2
    disc = b*b - 4 * a * c
    disc < 0 && return nothing
    sq = sqrt(disc)
    t1 = (-b - sq) / (2 * a)
    t2 = (-b + sq) / (2 * a)
    p1 = o + t1 * d
    p2 = o + t2 * d
    (Vector3{T}(p1[1], p1[2], p1[3]), Vector3{T}(p2[1], p2[2], p2[3]))
end

# ---------------------------------------------------------------------------
# Signed (oriented) angle between two vectors around an axis. The result lies
# in [0, 2π). Mirrors BALL's `getOrientedAngle`.
#
# Strips the component of v1 and v2 along `axis`, then measures the angle
# between the projected vectors in the plane perpendicular to `axis`.
# ---------------------------------------------------------------------------
function oriented_angle(v1::AbstractVector{T}, v2::AbstractVector{T}, axis::AbstractVector{T}) where {T<:Real}
    n = norm(axis)
    iszero(n) && return zero(T)
    a = axis / n
    p1 = v1 - dot(v1, a) * a
    p2 = v2 - dot(v2, a) * a
    n1 = norm(p1); n2 = norm(p2)
    (iszero(n1) || iszero(n2)) && return zero(T)
    cosθ = clamp(dot(p1, p2) / (n1 * n2), -one(T), one(T))
    sinθ = dot(cross(p1, p2), a) / (n1 * n2)
    θ = atan(sinθ, cosθ)
    θ < 0 ? θ + 2 * T(π) : θ
end

# ---------------------------------------------------------------------------
# Plane normal through three points (a, b, c). Returned vector is unit length;
# returns `nothing` if the points are collinear.
# ---------------------------------------------------------------------------
@inline function plane_normal(a::AbstractVector{T}, b::AbstractVector{T}, c::AbstractVector{T}) where {T<:Real}
    n = cross(b - a, c - a)
    ℓ = norm(n)
    iszero(ℓ) && return nothing
    Vector3{T}(n[1] / ℓ, n[2] / ℓ, n[3] / ℓ)
end

# ---------------------------------------------------------------------------
# Signed perpendicular distance from `p` to the plane through `a` with unit
# normal `n` (positive when `p` is on the side `n` points to).
# ---------------------------------------------------------------------------
@inline function plane_distance(p::AbstractVector{T}, a::AbstractVector{T}, n::AbstractVector{T}) where {T<:Real}
    dot(p - a, n)
end
