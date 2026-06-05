@testitem "icosphere" begin
    using LinearAlgebra
    for T in [Float32, Float64]
        # subdivisions = 0 → icosahedron
        m = icosphere(T, 0)
        @test nvertices(m) == 12
        @test ntriangles(m) == 20
        @test all(v -> isapprox(norm(v), 1; atol=sqrt(eps(T))), m.position)
        @test all(n -> isapprox(norm(n), 1; atol=sqrt(eps(T))), m.normal)

        # vertex count after n subdivisions: 10·4^n + 2
        for n in 1:3
            m = icosphere(T, n)
            @test nvertices(m) == 10 * 4^n + 2
            @test ntriangles(m) == 20 * 4^n
            @test all(v -> isapprox(norm(v), 1; atol=sqrt(eps(T))), m.position)
        end

        # area → 4π as subdivision increases (inscribed polyhedron, ~1/4ⁿ rate)
        @test surface_area(icosphere(T, 4)) ≈ 4*T(π) rtol = T(0.005)
    end
end

@testitem "surface_area" begin
    import GeometryBasics
    for T in [Float32, Float64]
        # Single triangle in the xy-plane with legs of length 3 and 4: area = 6.
        pos = [GeometryBasics.Point3{T}(0,0,0),
               GeometryBasics.Point3{T}(3,0,0),
               GeometryBasics.Point3{T}(0,4,0)]
        faces = [GeometryBasics.TriangleFace{Int}(1,2,3)]
        m = GeometryBasics.Mesh(pos, faces)
        @test surface_area(m) ≈ T(6)
    end
end

@testitem "Circle3" begin
    for T in [Float32, Float64]
        c = Circle3{T}(Vector3{T}(1, 0, 0), Vector3{T}(0, 0, 1), T(2))
        @test c isa Circle3{T}
        @test c.p == Vector3{T}(1, 0, 0)
        @test c.n == Vector3{T}(0, 0, 1)
        @test c.r == T(2)
    end
end

@testitem "read_msms" begin
    # Write a minimal MSMS-format pair to a temp dir so the test is
    # self-contained (no dependency on a local BALL checkout). MSMS
    # format: 3 header lines, then `x y z nx ny nz patch_id atom_id`
    # for .vert and `i j k face_type patch_id` for .face.
    mktempdir() do dir
        vpath = joinpath(dir, "test.vert")
        fpath = joinpath(dir, "test.face")
        open(vpath, "w") do io
            println(io, "# MSMS test fixture")
            println(io, "# generated for read_msms unit test")
            println(io, "3 2 0.0 0.0")
            println(io, "0.0 0.0 0.0  0.0 0.0 1.0  1 1 1")
            println(io, "1.0 0.0 0.0  0.0 0.0 1.0  1 1 1")
            println(io, "0.0 1.0 0.0  0.0 0.0 1.0  1 1 1")
        end
        open(fpath, "w") do io
            println(io, "# MSMS test fixture")
            println(io, "# generated for read_msms unit test")
            println(io, "1 1 0.0 0.0")
            println(io, "1 2 3 1 1")
        end
        for T in [Float32, Float64]
            m = read_msms(T, vpath, fpath)
            @test nvertices(m) == 3
            @test ntriangles(m) == 1
            @test m.position[1] ≈ T[0.0, 0.0, 0.0] atol=T(1e-3)
            @test m.normal[1]   ≈ T[0.0, 0.0, 1.0] atol=T(1e-3)
            @test Tuple(Int.(m.faces[1])) == (1, 2, 3)
        end
    end
end

@testitem "vdw_radius" begin
    for T in [Float32, Float64]
        @test vdw_radius(T, Elements.H) ≈ T(1.10)
        @test vdw_radius(T, Elements.C) ≈ T(1.70)
        @test vdw_radius(T, Elements.N) ≈ T(1.55)
        @test vdw_radius(T, Elements.O) ≈ T(1.52)
        @test vdw_radius(T, Elements.S) ≈ T(1.80)
        # fallback for elements without a tabulated value
        @test vdw_radius(T, Elements.Unknown) ≈ T(1.5)
        @test vdw_radius(T, Elements.C) isa T
    end
end

@testitem "assign_radii!" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)
        @test all(iszero, atoms(sys).radius)
        assign_radii!(sys)
        @test all(>(0), atoms(sys).radius)
        # carbons get 1.70, hydrogens 1.10, oxygens 1.52, nitrogens 1.55
        for (i, e) in enumerate(atoms(sys).element)
            @test atoms(sys).radius[i] ≈ vdw_radius(T, e)
        end

        # second call without overwrite is a no-op
        assign_radii!(sys)
        @test all(atoms(sys).radius .≈ vdw_radius.(T, atoms(sys).element))

        # overwrite=true replaces non-zero radii too
        atoms(sys).radius[1] = T(99)
        assign_radii!(sys; overwrite=true)
        @test atoms(sys).radius[1] ≈ vdw_radius(T, atoms(sys).element[1])
    end
end

@testitem "compute_numerical_sas: empty system" begin
    for T in [Float32, Float64]
        sys = System{T}()
        r = compute_numerical_sas(sys)
        @test r isa NumericalSASResult{T}
        @test r.total_area == zero(T)
        @test r.total_volume == zero(T)
        @test isempty(r.atom_areas)
        @test isempty(r.atom_volumes)
    end
end

@testitem "compute_numerical_sas: single atom matches 4π·R²" begin
    for T in [Float32, Float64]
        sys = System{T}()
        m = Molecule(sys)
        c = Chain(m)
        f = Fragment(c, 1; name="ALA")
        Atom(f, 1, Elements.C; r=Vector3{T}(0, 0, 0))
        assign_radii!(sys)

        probe = T(1.5)
        R = vdw_radius(T, Elements.C) + probe
        # Dense sampling — no occlusion → essentially exact.
        r = compute_numerical_sas(sys; probe_radius=probe, number_of_points=2000)
        @test r.total_area   ≈ 4*T(π)*R^2        rtol = T(1e-4)
        @test r.total_volume ≈ 4*T(π)*R^3/3      rtol = T(1e-4)
        @test length(r.atom_areas) == 1
        @test r.atom_areas[1]   ≈ r.total_area
        @test r.atom_volumes[1] ≈ r.total_volume
    end
end

@testitem "compute_numerical_sas: two distant atoms are additive" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        ch = Chain(mol)
        f = Fragment(ch, 1; name="ALA")
        Atom(f, 1, Elements.C; r=Vector3{T}(0, 0, 0))
        Atom(f, 2, Elements.C; r=Vector3{T}(20, 0, 0))
        assign_radii!(sys)

        probe = T(1.5)
        R = vdw_radius(T, Elements.C) + probe
        r = compute_numerical_sas(sys; probe_radius=probe, number_of_points=2000)
        @test r.total_area   ≈ 2 * 4*T(π)*R^2   rtol = T(1e-4)
        @test r.total_volume ≈ 2 * 4*T(π)*R^3/3 rtol = T(1e-4)
        @test r.atom_areas[1]   ≈ r.atom_areas[2]   rtol = T(1e-4)
        @test r.atom_volumes[1] ≈ r.atom_volumes[2] rtol = T(1e-4)
    end
end

@testitem "compute_numerical_sas: overlapping atoms reduce area" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        ch = Chain(mol)
        f = Fragment(ch, 1; name="ALA")
        Atom(f, 1, Elements.C; r=Vector3{T}(0, 0, 0))
        Atom(f, 2, Elements.C; r=Vector3{T}(1.5, 0, 0))
        assign_radii!(sys)

        probe = T(1.5)
        R = vdw_radius(T, Elements.C) + probe
        r = compute_numerical_sas(sys; probe_radius=probe, number_of_points=2000)
        # Two heavily overlapping spheres → area must be strictly less than
        # the additive case but more than a single sphere.
        @test r.total_area < 2 * 4*T(π)*R^2
        @test r.total_area > 4*T(π)*R^2
        # By symmetry, per-atom areas are equal.
        @test r.atom_areas[1] ≈ r.atom_areas[2] rtol = T(1e-3)
    end
end

@testitem "compute_numerical_sas: skips zero-radius atoms" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        ch = Chain(mol)
        f = Fragment(ch, 1; name="ALA")
        Atom(f, 1, Elements.C; r=Vector3{T}(0, 0, 0))
        Atom(f, 2, Elements.H; r=Vector3{T}(0, 5, 0))
        # only assign the carbon's radius — hydrogen stays at 0
        atoms(sys).radius[1] = vdw_radius(T, Elements.C)

        r = compute_numerical_sas(sys; probe_radius=T(1.5), number_of_points=500)
        @test r.atom_areas[2] == 0
        @test r.atom_volumes[2] == 0
        # carbon's area is the isolated-sphere area
        R = vdw_radius(T, Elements.C) + T(1.5)
        @test r.total_area ≈ 4*T(π)*R^2 rtol = T(5e-3)
    end
end

@testitem "compute_numerical_sas: surface point cloud" begin
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys)
        ch = Chain(mol)
        f = Fragment(ch, 1; name="ALA")
        Atom(f, 1, Elements.C; r=Vector3{T}(0, 0, 0))
        assign_radii!(sys)

        probe = T(1.5)
        R = vdw_radius(T, Elements.C) + probe
        r = compute_numerical_sas(sys; probe_radius=probe, number_of_points=200,
                                   compute_surface=true, compute_surface_per_atom=true)
        # all points on the sphere surface
        @test all(v -> isapprox(sqrt(sum(abs2, v)), R; atol=sqrt(eps(T))*R), r.surface_vertices)
        # the sum of normal lengths equals total area
        @test sum(v -> sqrt(sum(abs2, v)), r.surface_normals) ≈ r.total_area rtol=T(1e-4)
        @test length(r.atom_surface_vertices) == 1
        @test length(r.atom_surface_vertices[1]) == length(r.surface_vertices)
    end
end

@testitem "compute_numerical_sas: bpti smoke test" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        r = compute_numerical_sas(sys; probe_radius=T(1.5), number_of_points=400)
        # BPTI's SAS area is published in the 3500-4500 Å² range with these
        # radii and probe; we sample relatively coarse, so allow a wide window.
        @test 3000 < r.total_area  < 5000
        @test 8000 < r.total_volume < 14000
        @test length(r.atom_areas) == natoms(sys)
        @test sum(r.atom_areas)   ≈ r.total_area   rtol = T(1e-5)
        @test sum(r.atom_volumes) ≈ r.total_volume rtol = T(1e-5)
    end
end

@testitem "compute_reduced_surface: empty / trivial systems" begin
    for T in [Float32, Float64]
        # Empty
        let rs = compute_reduced_surface(Sphere{T}[]; probe_radius=T(1.0))
            @test rs isa ReducedSurface{T}
            @test length(rs.atoms) == 0
            @test isempty(rs.vertices) && isempty(rs.edges) && isempty(rs.faces)
        end
        # Single atom — no face, but the atom is registered as an isolated vertex
        let rs = compute_reduced_surface(
                [Sphere{T}(BiochemicalAlgorithms.Point3{T}(0,0,0), T(1.5))];
                probe_radius=T(1.0))
            @test isempty(rs.faces)
            @test BiochemicalAlgorithms.nvertices(rs) == 1
        end
    end
end

@testitem "compute_reduced_surface: equilateral triangle" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        # Three identical atoms in the xy-plane at the vertices of an
        # equilateral triangle far enough apart that the probe rests
        # tangentially on all three without occlusion.
        a = T(3)
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a/2,       a*sqrt(T(3))/2, 0), T(1.5)),
        ]
        rs = compute_reduced_surface(spheres; probe_radius=T(1.0))
        # Two faces (probe above and below the plane), 3 vertices (one per
        # atom), 3 edges (one per atom-pair, connecting the two faces).
        @test BiochemicalAlgorithms.nfaces(rs) == 2
        @test BiochemicalAlgorithms.nvertices(rs) == 3
        @test BiochemicalAlgorithms.nedges(rs) == 3
        # 3F = 2E for a closed manifold
        @test 3 * BiochemicalAlgorithms.nfaces(rs) == 2 * BiochemicalAlgorithms.nedges(rs)
        # Both faces should have all three edges populated and reference
        # the same edge set.
        for f in rs.faces
            f.v1 == 0 && continue   # sentinel = deleted face
            @test length(f.edges) >= 3
            @test all(>(0), f.edges)
        end
        # All edges have both adjacent faces set
        for e in rs.edges
            e.v1 == 0 && continue   # sentinel = deleted edge
            @test e.f1 != 0 && e.f2 != 0
        end
        # The two probe centers should be mirror images across the xy plane.
        zs = sort!([f.center[3] for f in rs.faces])
        @test zs[1] ≈ -zs[2] atol=sqrt(eps(T))
        # Normals point outward, away from the plane → opposite z signs.
        ns = sort!([f.normal[3] for f in rs.faces])
        @test ns[1] < 0 < ns[2]
    end
end

@testitem "compute_reduced_surface: two-atom system has a free edge" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        # Probe touching two close-enough atoms can rotate freely → no face
        # is born, but a single free edge with full 2π rolling is recorded.
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,0,0), T(1.5)),
            Sphere{T}(Point3{T}(3,0,0), T(1.5)),
        ]
        rs = compute_reduced_surface(spheres; probe_radius=T(1.0))
        @test isempty(rs.faces)
        @test BiochemicalAlgorithms.nedges(rs) == 1
        e = only(rs.edges)
        @test e.f1 == 0 && e.f2 == 0
        @test e.angle ≈ 2 * T(π)
        @test BiochemicalAlgorithms.nvertices(rs) == 2
    end
end

@testitem "compute_reduced_surface: isolated single atom yields one vertex" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        rs = compute_reduced_surface(
            [Sphere{T}(Point3{T}(0,0,0), T(1.5))]; probe_radius=T(1.0))
        @test BiochemicalAlgorithms.nvertices(rs) == 1
        @test isempty(rs.edges)
        @test isempty(rs.faces)
    end
end

@testitem "compute_reduced_surface: distant atoms yield isolated vertices" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        # Two atoms 20 Å apart — probe can't reach both at once.
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,0,0), T(1.5)),
            Sphere{T}(Point3{T}(20,0,0), T(1.5)),
        ]
        rs = compute_reduced_surface(spheres; probe_radius=T(1.0))
        @test isempty(rs.faces)
        @test isempty(rs.edges)
        @test BiochemicalAlgorithms.nvertices(rs) == 2
    end
end

@testitem "compute_reduced_surface: bpti smoke test" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        rs = compute_reduced_surface(sys; probe_radius=T(1.5))
        @test rs isa ReducedSurface{T}
        @test BiochemicalAlgorithms.nfaces(rs) > 0
        @test BiochemicalAlgorithms.nvertices(rs) > 0
        @test BiochemicalAlgorithms.nedges(rs) > 0
        # Every face references exactly 3 distinct vertices on 3 distinct atoms.
        for f in rs.faces
            atoms_set = Set((rs.vertices[f.v1].atom,
                             rs.vertices[f.v2].atom,
                             rs.vertices[f.v3].atom))
            @test length(atoms_set) == 3
        end
        # Probe radius and bounding box are populated.
        @test rs.probe_radius == T(1.5)
        @test rs.r_max > 0
    end
end

@testitem "compute_sas: SAS is the dual of RS" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        a = T(3)
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a/2,       a*sqrt(T(3))/2, 0), T(1.5)),
        ]
        rs  = compute_reduced_surface(spheres; probe_radius=T(1.0))
        sas = compute_sas(rs)
        @test sas isa SolventAccessibleSurface{T}
        # Duality: V_sas = F_rs, E_sas = E_rs, F_sas = V_rs
        @test BiochemicalAlgorithms.nvertices(sas) == BiochemicalAlgorithms.nfaces(rs)
        @test BiochemicalAlgorithms.nedges(sas)    == BiochemicalAlgorithms.nedges(rs)
        @test BiochemicalAlgorithms.nfaces(sas)    == BiochemicalAlgorithms.nvertices(rs)
        # Each SAS vertex sits at the center of the RS probe sphere.
        # SAS arrays are index-aligned with RS arrays and carry sentinels
        # for entries whose RS counterpart was removed during construction
        # (rs.faces[j].v1 == 0 / rs.vertices[j].atom == 0); skip those.
        for (j, v) in enumerate(sas.vertices)
            rs.faces[j].v1 == 0 && continue
            @test v.point ≈ rs.faces[j].center
        end
        # Each SAS face's sphere is the inflated atom (radius + probe_radius).
        for (j, f) in enumerate(sas.faces)
            rs.vertices[j].atom == 0 && continue
            atom = rs.atoms[rs.vertices[j].atom]
            @test f.sphere.r ≈ atom.r + rs.probe_radius
            @test f.sphere.center ≈ atom.center
        end
    end
end

@testitem "compute_sas: bpti smoke test" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        sas = compute_sas(sys; probe_radius=T(1.5))
        @test BiochemicalAlgorithms.nvertices(sas) > 0
        @test BiochemicalAlgorithms.nedges(sas)    > 0
        @test BiochemicalAlgorithms.nfaces(sas)    > 0
        @test BiochemicalAlgorithms.nvertices(sas) == BiochemicalAlgorithms.nfaces(sas.reduced_surface)
    end
end

@testitem "compute_ses: SES face counts match RS" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        a = T(3)
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a/2,       a*sqrt(T(3))/2, 0), T(1.5)),
        ]
        rs  = compute_reduced_surface(spheres; probe_radius=T(1.0))
        ses = compute_ses(rs)
        @test ses isa SolventExcludedSurface{T}
        # Face-count duality
        @test BiochemicalAlgorithms.nspheric_faces(ses) == BiochemicalAlgorithms.nfaces(rs)
        @test BiochemicalAlgorithms.ntoric_faces(ses)   == BiochemicalAlgorithms.nedges(rs)
        @test BiochemicalAlgorithms.ncontact_faces(ses) == BiochemicalAlgorithms.nvertices(rs)
        # Three contact points per spheric face
        @test BiochemicalAlgorithms.nvertices(ses) == 3 * BiochemicalAlgorithms.nfaces(rs)
        # All spheric face vertices lie on the probe sphere of their RS face.
        for spf in ses.spheric_faces
            for vi in spf.vertices
                p = ses.vertices[vi].point
                @test isapprox(sqrt(sum(abs2, p - spf.sphere.center)),
                               spf.sphere.r; atol = sqrt(eps(T)) * 100)
            end
        end
        # Contact points lie on their atom's surface
        for v in ses.vertices
            atom = rs.atoms[v.atom]
            @test isapprox(sqrt(sum(abs2, v.point - atom.center)),
                           atom.r; atol = sqrt(eps(T)) * 100)
        end
    end
end

@testitem "compute_ses: bpti smoke test" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        ses = compute_ses(sys; probe_radius=T(1.5))
        @test BiochemicalAlgorithms.nspheric_faces(ses) > 0
        @test BiochemicalAlgorithms.ntoric_faces(ses)   > 0
        @test BiochemicalAlgorithms.ncontact_faces(ses) > 0
        # Each SES vertex points back to either an RS face (regular contact
        # point) or an RS edge (singular toric intersection point).
        for v in ses.vertices
            face_ok = 1 <= v.rs_face <= length(ses.reduced_surface.faces)
            edge_ok = 1 <= v.rs_edge <= length(ses.reduced_surface.edges)
            @test face_ok || edge_ok
        end
    end
end

@testitem "triangulate_sas: produces a valid mesh" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        a = T(3)
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a/2,       a*sqrt(T(3))/2, 0), T(1.5)),
        ]
        sas = compute_sas(compute_reduced_surface(spheres; probe_radius=T(1.0)))
        m = triangulate_sas(sas)
        @test nvertices(m) > 0
        @test ntriangles(m) > 0
        # All vertices on one of the inflated atom spheres (radius 2.5).
        for v in m.position
            ok = any(sas.faces) do f
                d = sqrt(sum(abs2, Vector{T}(v) - Vector{T}(f.sphere.center)))
                isapprox(d, f.sphere.r; atol = T(1e-4))
            end
            @test ok
        end
        # Higher density gives more triangles. Surface area is non-
        # monotonic: BALL's SAS triangulator drops every triangle whose
        # corners are all behind a neighbour cut plane, so denser sampling
        # increases the number of clipped triangles at boundaries — net
        # area can decrease as density rises (see memory:
        # `surfaces_knowledge` → "BALL's SAS is leaky by design").
        m1 = triangulate_sas(sas; density=T(1))
        m2 = triangulate_sas(sas; density=T(10))
        @test ntriangles(m2) > ntriangles(m1)
    end
end

@testitem "triangulate_sas: bpti smoke test" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        m = triangulate_sas(sys; probe_radius=T(1.5), density=T(2))
        @test nvertices(m) > 1000
        @test ntriangles(m) > 1000
        @test surface_area(m) > 1000
    end
end

@testitem "triangulate_ses: produces a valid mesh" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        a = T(3)
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a,         0,             0), T(1.5)),
            Sphere{T}(Point3{T}(a/2,       a*sqrt(T(3))/2, 0), T(1.5)),
        ]
        ses = compute_ses(compute_reduced_surface(spheres; probe_radius=T(1.0)))
        m = triangulate_ses(ses)
        @test nvertices(m) > 0
        @test ntriangles(m) > 0
        # Higher density: more triangles.
        m1 = triangulate_ses(ses; density=T(1))
        m2 = triangulate_ses(ses; density=T(10))
        @test ntriangles(m2) > ntriangles(m1)
    end
end

@testitem "triangulate_ses: bpti smoke test" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        m = triangulate_ses(sys; probe_radius=T(1.5), density=T(2))
        @test nvertices(m) > 1000
        @test ntriangles(m) > 1000
        @test surface_area(m) > 1000
    end
end

@testitem "triangulate: dispatches on surface type" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        a = T(3)
        spheres = Sphere{T}[
            Sphere{T}(Point3{T}(0,0,0), T(1.5)),
            Sphere{T}(Point3{T}(a,0,0), T(1.5)),
            Sphere{T}(Point3{T}(a/2, a*sqrt(T(3))/2, 0), T(1.5)),
        ]
        rs  = compute_reduced_surface(spheres; probe_radius=T(1.0))
        sas = compute_sas(rs)
        ses = compute_ses(rs)
        m_sas = triangulate(sas)
        m_ses = triangulate(ses)
        @test ntriangles(m_sas) == ntriangles(triangulate_sas(sas))
        @test ntriangles(m_ses) == ntriangles(triangulate_ses(ses))
    end
end

@testitem "sas_area / sas_volume: single atom is exact" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        sys = System{T}()
        m = Molecule(sys); c = Chain(m); f = Fragment(c, 1; name="X")
        Atom(f, 1, Elements.C; r=Vector3{T}(0,0,0))
        assign_radii!(sys)
        R = vdw_radius(T, Elements.C) + T(1.5)
        @test sas_area(sys; number_of_points=2000)   ≈ 4*T(π)*R^2  rtol=T(1e-3)
        @test sas_volume(sys; number_of_points=2000) ≈ 4*T(π)*R^3/3 rtol=T(1e-3)
    end
end

@testitem "ses_area: single atom is exact" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        sys = System{T}()
        m = Molecule(sys); c = Chain(m); f = Fragment(c, 1; name="X")
        Atom(f, 1, Elements.C; r=Vector3{T}(0,0,0))
        assign_radii!(sys)
        r = vdw_radius(T, Elements.C)
        # Isolated atom: SES = bare atom area, no toric/spheric contributions.
        @test ses_area(sys; number_of_points=2000) ≈ 4*T(π)*r^2 rtol=T(1e-3)
    end
end

@testitem "ses_area: matches numerical SAS scale for BPTI" begin
    # The analytical SES sum (contact + toric + spheric) should be in the
    # same order of magnitude as published BPTI values (≈ 3500-9000 Å²
    # depending on the singularity-handling and probe-radius choice).
    # We just sanity-check that it's positive and within a wide window.
    T = Float64
    sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
    assign_radii!(sys)
    A_ses = ses_area(sys; probe_radius=T(1.5), number_of_points=400)
    @test 2000 < A_ses < 12000
end

@testitem "check_rs / check_ses: topology of a clean BPTI run is healthy" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        rs = compute_reduced_surface(sys; probe_radius=T(1.5))
        ses = compute_ses(rs)
        @test isempty(check_rs(rs))
        @test isempty(check_ses(ses))
    end
end

@testitem "clean_ses! / split_spheric_faces!: idempotent on healthy SES" begin
    using BiochemicalAlgorithms: Point3
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        ses = compute_ses(sys; probe_radius=T(1.5))
        # default thresholds: nothing to drop
        @test clean_ses!(ses) == 0
        # No splitting needed for a typical molecule
        @test split_spheric_faces!(ses) == 0
        @test isempty(check_ses(ses))
        # Aggressive threshold drops some toric faces
        ses2 = compute_ses(sys; probe_radius=T(1.5))
        dropped = clean_ses!(ses2; min_angle=T(0.5))
        @test dropped > 0
        @test isempty(check_ses(ses2))
    end
end

@testitem "resolve_probe_intersections!: detects overlapping probes" begin
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/bpti.pdb"), T)
        assign_radii!(sys)
        ses = compute_ses(sys; probe_radius=T(1.5))
        n_before = length(ses.vertices)
        n_pairs = resolve_probe_intersections!(ses)
        @test n_pairs > 0
        # Each handled pair adds two new SES vertices
        @test length(ses.vertices) == n_before + 2 * n_pairs
        @test isempty(check_ses(ses))
    end
end

@testitem "compute_reduced_surface: detects spindle-torus singularities" begin
    T = Float64
    sys = System{T}()
    mol = Molecule(sys); ch = Chain(mol); f = Fragment(ch, 1; name="X")
    # Square of four C atoms; opposite-diagonal pair has spindle torus.
    Atom(f, 1, Elements.C; r=Vector3{T}(0,    0,    0))
    Atom(f, 2, Elements.C; r=Vector3{T}(5.8,  0,    0))
    Atom(f, 3, Elements.C; r=Vector3{T}(2.9,  2.9,  0))
    Atom(f, 4, Elements.C; r=Vector3{T}(2.9, -2.9,  0))
    for i in 1:4
        atoms(sys).radius[i] = T(1.7)
    end
    rs = compute_reduced_surface(sys; probe_radius=T(1.5))
    @test any(e -> e.singular, rs.edges)
    # For singular edges the intersection points lie on the probe sphere of
    # the corresponding RS face.
    for e in rs.edges
        e.singular || continue
        e.f1 != 0 || continue
        probe = rs.faces[e.f1].center
        d1 = e.intersection1 - probe
        d2 = e.intersection2 - probe
        @test isapprox(sqrt(sum(abs2, d1)), rs.probe_radius; atol=T(1e-4))
        @test isapprox(sqrt(sum(abs2, d2)), rs.probe_radius; atol=T(1e-4))
    end
    ses = compute_ses(rs)
    @test any(f -> f.type == SESFaceType.ToricSingular, ses.toric_faces)
end
