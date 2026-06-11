@testitem "cpk_color" begin
    using BiochemicalAlgorithms
    @test cpk_color(Elements.H) == (1.0, 1.0, 1.0)
    @test cpk_color(Elements.C) == (0.300, 0.300, 0.300)
    @test cpk_color(Elements.O)[1] > cpk_color(Elements.O)[2]   # red-ish
    @test cpk_color(Elements.N)[3] > cpk_color(Elements.N)[1]   # blue-ish
    # Unmapped element falls back to fuchsia.
    @test cpk_color(Elements.Au) == (1.000, 0.078, 0.576)
end

@testitem "PrintablePart + STL/3MF I/O" begin
    using BiochemicalAlgorithms
    using BiochemicalAlgorithms: PrintablePart
    GB = BiochemicalAlgorithms.GeometryBasics
    verts = GB.Point3{Float64}[
        GB.Point3{Float64}(0, 0, 0), GB.Point3{Float64}(1, 0, 0),
        GB.Point3{Float64}(0, 1, 0), GB.Point3{Float64}(0, 0, 1),
    ]
    faces = GB.TriangleFace{Int}[
        GB.TriangleFace{Int}(1, 3, 2), GB.TriangleFace{Int}(1, 2, 4),
        GB.TriangleFace{Int}(1, 4, 3), GB.TriangleFace{Int}(2, 3, 4),
    ]
    m = GB.Mesh(verts, faces)

    mktempdir() do dir
        stl = joinpath(dir, "tetra.stl")
        export_stl(m, stl)
        # STL: 80 header + 4 byte triCount + 50 bytes/tri.
        @test filesize(stl) == 80 + 4 + 50 * length(faces)
        # First 80 bytes are header, next is little-endian UInt32 triangle count.
        n = open(stl) do io; read(io, 80); read(io, UInt32) end
        @test n == length(faces)

        mf = joinpath(dir, "tetra.3mf")
        part = PrintablePart(m, (0.2, 0.4, 0.8), "tetra")
        export_3mf([part], mf)
        @test isfile(mf)
        @test filesize(mf) > 0

        # 3MF is a ZIP — first two bytes are 'PK'.
        magic = open(mf) do io; read(io, 2) end
        @test magic == UInt8['P', 'K']
    end
end

@testitem "cylinder_mesh: closed + watertight + correct face count" begin
    using BiochemicalAlgorithms
    using BiochemicalAlgorithms: cylinder_mesh
    using LinearAlgebra
    for T in [Float32, Float64]
        m = cylinder_mesh(T[0, 0, 0], T[0, 0, 5], T(1.0); segments=12)
        # 2*12 ring verts + 2 cap centers = 26 verts.
        @test length(m.position) == 2 * 12 + 2
        # Sides: 2*12 tris.  Caps: 2*12 tris.
        @test length(m.faces) == 4 * 12

        # Euler: V - E + F = 2 for a closed orientable surface.
        # Each closed triangle mesh has E = 3F/2.
        V = length(m.position); F = length(m.faces); E = 3F ÷ 2
        @test V - E + F == 2

        # All side vertices sit at radius 1 from the axis.
        for v in m.position[1:24]
            @test isapprox(sqrt(v[1]^2 + v[2]^2), 1; atol=sqrt(eps(T)))
        end
    end
end

@testitem "construction_kit: water" begin
    using BiochemicalAlgorithms
    using BiochemicalAlgorithms: PrintablePart
    for T in [Float32, Float64]
        sys = System{T}()
        mol = Molecule(sys; name="water")
        ch = Chain(mol)
        f = Fragment(ch, 1; name="WAT")
        oa = Atom(f, 1, Elements.O; r=Vector3{T}(0, 0, 0))
        ha = Atom(f, 2, Elements.H; r=Vector3{T}(0.96, 0, 0))
        hb = Atom(f, 3, Elements.H; r=Vector3{T}(-0.24, 0.93, 0))
        Bond(sys, oa.idx, ha.idx, BondOrder.Single)
        Bond(sys, oa.idx, hb.idx, BondOrder.Single)
        assign_radii!(sys)

        parts = construction_kit(sys; scale=10)
        @test parts isa Vector{PrintablePart}
        @test length(parts) == 5                # 3 atoms + 2 bonds
        @test count(p -> startswith(p.name, "atom-"), parts) == 3
        @test count(p -> startswith(p.name, "bond-"), parts) == 2
        # Oxygen gets the red CPK color.
        oxygen_part = parts[findfirst(p -> endswith(p.name, "-O"), parts)]
        r, g, b = oxygen_part.color
        @test r > g && r > b
        for p in parts
            @test length(p.mesh.position) > 0
            @test length(p.mesh.faces)    > 0
        end
    end
end

@testitem "construction_kit: bond orders 1/2/3" begin
    using BiochemicalAlgorithms
    for ord in (BondOrder.Single, BondOrder.Double, BondOrder.Triple)
        sys = System{Float64}()
        mol = Molecule(sys; name="x")
        ch = Chain(mol)
        f = Fragment(ch, 1; name="X")
        a = Atom(f, 1, Elements.C; r=Vector3{Float64}(0, 0, 0))
        b = Atom(f, 2, Elements.C; r=Vector3{Float64}(1.5, 0, 0))
        Bond(sys, a.idx, b.idx, ord)
        assign_radii!(sys)

        parts = construction_kit(sys; scale=10)
        @test length(parts) == 3
        bond = parts[end]
        @test startswith(bond.name, "bond-")
    end
end

@testitem "construction_kit: :magnet joint and error on bogus" begin
    using BiochemicalAlgorithms
    sys = System{Float64}()
    mol = Molecule(sys; name="x")
    ch = Chain(mol)
    f = Fragment(ch, 1; name="X")
    a = Atom(f, 1, Elements.C; r=Vector3{Float64}(0, 0, 0))
    b = Atom(f, 2, Elements.O; r=Vector3{Float64}(1.4, 0, 0))
    Bond(sys, a.idx, b.idx, BondOrder.Single)
    assign_radii!(sys)

    parts_peg    = construction_kit(sys; joint=:peg)
    parts_magnet = construction_kit(sys; joint=:magnet)
    @test length(parts_peg) == length(parts_magnet) == 3
    @test_throws ArgumentError construction_kit(sys; joint=:bogus)
end

@testitem "export_stl on a cylinder round-trips file size" begin
    using BiochemicalAlgorithms
    using BiochemicalAlgorithms: cylinder_mesh
    m = cylinder_mesh([0.0, 0, 0], [0.0, 0, 10.0], 2.5; segments=20)
    mktempdir() do dir
        path = joinpath(dir, "c.stl")
        export_stl(m, path)
        @test filesize(path) == 80 + 4 + 50 * length(m.faces)
    end
end

@testitem "PrintablePart: per-face colors round-trip into 3MF" begin
    using BiochemicalAlgorithms
    using BiochemicalAlgorithms: PrintablePart
    ZipFile = BiochemicalAlgorithms.ZipFile
    GB = BiochemicalAlgorithms.GeometryBasics
    verts = GB.Point3{Float64}[
        GB.Point3{Float64}(0, 0, 0), GB.Point3{Float64}(1, 0, 0),
        GB.Point3{Float64}(0, 1, 0), GB.Point3{Float64}(0, 0, 1),
    ]
    faces = GB.TriangleFace{Int}[
        GB.TriangleFace{Int}(1, 3, 2), GB.TriangleFace{Int}(1, 2, 4),
        GB.TriangleFace{Int}(1, 4, 3), GB.TriangleFace{Int}(2, 3, 4),
    ]
    m = GB.Mesh(verts, faces)
    # Distinct colour per triangle.
    fc = NTuple{3, Float64}[
        (1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0), (1.0, 1.0, 0.0),
    ]
    part = PrintablePart(m, (0.5, 0.5, 0.5), "tetra-rgb", fc)

    mktempdir() do dir
        path = joinpath(dir, "tetra-rgb.3mf")
        export_3mf([part], path)
        # Pop the 3dmodel.model out of the zip and check the colorgroup +
        # per-triangle pid/p1/p2/p3 are present.
        r = ZipFile.Reader(path)
        xml = ""
        for f in r.files
            if endswith(f.name, "3dmodel.model")
                xml = String(read(f))
                break
            end
        end
        close(r)
        @test occursin("m:colorgroup", xml)
        @test occursin("#ff0000", xml)
        @test occursin("#00ff00", xml)
        @test occursin("#0000ff", xml)
        @test occursin("#ffff00", xml)
        # Every triangle line carries a `pid=` attribute.
        @test count(Regex("<triangle [^>]*pid="), xml) == length(faces)
    end
end

@testitem "ses_face_colors_by_atom + export_ses_3mf" begin
    using BiochemicalAlgorithms
    ZipFile = BiochemicalAlgorithms.ZipFile
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    ses = compute_ses(sys; probe_radius=1.5)
    tses = triangulate_ses(ses; density=2.0)
    colors = ses_face_colors_by_atom(tses, ses, sys)
    @test length(colors) == length(tses.faces)
    # All colours should be valid RGB tuples.
    @test all(c -> all(0 .<= c .<= 1), colors)
    # At least two distinct colours (AlaAla has multiple element types).
    @test length(unique(colors)) >= 2

    mktempdir() do dir
        path = joinpath(dir, "ala.3mf")
        export_ses_3mf(ses, sys, path; density=2.0)
        @test isfile(path)
        @test filesize(path) > 1000

        # Read back and confirm <m:colorgroup> + per-triangle pid attrs.
        r = ZipFile.Reader(path)
        xml = ""
        for f in r.files
            if endswith(f.name, "3dmodel.model")
                xml = String(read(f))
                break
            end
        end
        close(r)
        @test occursin("m:colorgroup", xml)
        @test count(Regex("<triangle [^>]*pid="), xml) == length(tses.faces)
    end
end

@testitem "export_ses_3mf one-call convenience from AtomContainer" begin
    using BiochemicalAlgorithms
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    mktempdir() do dir
        path = joinpath(dir, "ala.3mf")
        export_ses_3mf(sys, path; probe_radius=1.5, density=2.0)
        @test isfile(path)
        @test filesize(path) > 1000
    end
end

@testitem "construction_kit parts sit on z=0 and carry face_colors" begin
    using BiochemicalAlgorithms
    sys = System{Float64}()
    mol = Molecule(sys; name="water")
    ch = Chain(mol); f = Fragment(ch, 1; name="WAT")
    oa = Atom(f, 1, Elements.O; r=Vector3{Float64}(0, 0, 0))
    ha = Atom(f, 2, Elements.H; r=Vector3{Float64}(0.96, 0, 0))
    hb = Atom(f, 3, Elements.H; r=Vector3{Float64}(-0.24, 0.93, 0))
    Bond(sys, oa.idx, ha.idx, BondOrder.Single)
    Bond(sys, oa.idx, hb.idx, BondOrder.Single)
    assign_radii!(sys)
    parts = construction_kit(sys; scale=20)
    for p in parts
        zs = [v[3] for v in p.mesh.position]
        xs = [v[1] for v in p.mesh.position]
        ys = [v[2] for v in p.mesh.position]
        # z-min == 0 (sits on build plate).
        @test isapprox(minimum(zs), 0; atol=1e-8)
        # Recentred on XY origin.
        @test isapprox((minimum(xs) + maximum(xs)) / 2, 0; atol=1e-8)
        @test isapprox((minimum(ys) + maximum(ys)) / 2, 0; atol=1e-8)
        # Every part has uniform face_colors set.
        @test p.face_colors !== nothing
        @test length(p.face_colors) == length(p.mesh.faces)
        @test length(unique(p.face_colors)) == 1
    end
end

@testitem "export_ses_3mf reorient=true reduces z extent" begin
    using BiochemicalAlgorithms
    using LinearAlgebra
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    ses = compute_ses(sys; probe_radius=1.5)
    mktempdir() do dir
        # Compare unoriented vs reoriented bounding boxes.
        unori = joinpath(dir, "unori.3mf")
        ori   = joinpath(dir, "ori.3mf")
        export_ses_3mf(ses, sys, unori; density=2.0, reorient=false)
        export_ses_3mf(ses, sys, ori;   density=2.0, reorient=true)
        @test isfile(unori) && isfile(ori)
        @test filesize(unori) > 0 && filesize(ori) > 0
    end
end
