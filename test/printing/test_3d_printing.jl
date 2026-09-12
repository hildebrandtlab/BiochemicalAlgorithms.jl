@testitem "cpk_color" begin
    @test cpk_color(Elements.H) == (1.0, 1.0, 1.0)
    @test cpk_color(Elements.C) == (0.300, 0.300, 0.300)
    @test cpk_color(Elements.O)[1] > cpk_color(Elements.O)[2]   # red-ish
    @test cpk_color(Elements.N)[3] > cpk_color(Elements.N)[1]   # blue-ish
    # Unmapped element falls back to fuchsia.
    @test cpk_color(Elements.Au) == (1.000, 0.078, 0.576)
end

@testitem "PrintablePart + STL/3MF I/O" begin
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

@testitem "construction_kit: bond orders 1/2/3 emit order-many bond parts" begin
    # Single bond → 1 bond piece; double → 2 curved cylinders; triple → 3.
    for (ord, n_bond_parts) in (
        (BondOrder.Single, 1),
        (BondOrder.Double, 2),
        (BondOrder.Triple, 3),
    )
        sys = System{Float64}()
        mol = Molecule(sys; name="x")
        ch = Chain(mol); f = Fragment(ch, 1; name="X")
        a = Atom(f, 1, Elements.C; r=Vector3{Float64}(0, 0, 0))
        b = Atom(f, 2, Elements.C; r=Vector3{Float64}(1.5, 0, 0))
        Bond(sys, a.idx, b.idx, ord)
        assign_radii!(sys)

        parts = construction_kit(sys; scale=10)
        bond_parts = filter(p -> startswith(p.name, "bond-"), parts)
        @test length(parts) == 2 + n_bond_parts
        @test length(bond_parts) == n_bond_parts
    end
end

@testitem "construction_kit: curved double-bond cylinders are watertight" begin
    function audit(mesh)
        V = length(mesh.position); F = length(mesh.faces)
        ec = Dict{Tuple{Int,Int}, Int}()
        for fa in mesh.faces
            for (a, b) in ((Int(fa[1]),Int(fa[2])),(Int(fa[2]),Int(fa[3])),(Int(fa[3]),Int(fa[1])))
                k = a < b ? (a,b) : (b,a); ec[k] = get(ec, k, 0) + 1
            end
        end
        E = length(ec); bd = count(==(1), values(ec)); nm = count(>(2), values(ec))
        (V=V, E=E, F=F, χ=V-E+F, boundary=bd, non_manifold=nm)
    end

    # Ethylene C=C: classic curved-double-bond test case
    sys = System{Float64}()
    mol = Molecule(sys; name="ethylene")
    ch = Chain(mol); f = Fragment(ch, 1; name="ETH")
    c1 = Atom(f, 1, Elements.C; r=Vector3{Float64}( 0.67, 0, 0))
    c2 = Atom(f, 2, Elements.C; r=Vector3{Float64}(-0.67, 0, 0))
    h1 = Atom(f, 3, Elements.H; r=Vector3{Float64}( 1.24,  0.92, 0))
    h2 = Atom(f, 4, Elements.H; r=Vector3{Float64}( 1.24, -0.92, 0))
    h3 = Atom(f, 5, Elements.H; r=Vector3{Float64}(-1.24,  0.92, 0))
    h4 = Atom(f, 6, Elements.H; r=Vector3{Float64}(-1.24, -0.92, 0))
    Bond(sys, c1.idx, c2.idx, BondOrder.Double)
    for h in (h1, h2); Bond(sys, c1.idx, h.idx, BondOrder.Single); end
    for h in (h3, h4); Bond(sys, c2.idx, h.idx, BondOrder.Single); end
    assign_radii!(sys)
    parts = construction_kit(sys; scale=20)
    # 2 C + 4 H + (2 double-bond cylinders) + (4 C-H bonds) = 12 parts
    @test length(parts) == 12
    # The two double-bond cylinders are named bond-1-cyl-1 and bond-1-cyl-2
    cyls = filter(p -> startswith(p.name, "bond-1-cyl-"), parts)
    @test length(cyls) == 2
    for p in cyls
        s = audit(p.mesh)
        @test s.χ == 2
        @test s.boundary == 0
        @test s.non_manifold == 0
    end
end

@testitem "construction_kit: :magnet joint and error on bogus" begin
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
    using BiochemicalAlgorithms: cylinder_mesh
    m = cylinder_mesh([0.0, 0, 0], [0.0, 0, 10.0], 2.5; segments=20)
    mktempdir() do dir
        path = joinpath(dir, "c.stl")
        export_stl(m, path)
        @test filesize(path) == 80 + 4 + 50 * length(m.faces)
    end
end

@testitem "PrintablePart: per-face colors round-trip into 3MF" begin
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
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    mktempdir() do dir
        path = joinpath(dir, "ala.3mf")
        export_ses_3mf(sys, path; probe_radius=1.5, density=2.0)
        @test isfile(path)
        @test filesize(path) > 1000
    end
end

@testitem "export_ses_3mf max_size_mm scales mesh to target build-volume" begin
    ZipFile = BiochemicalAlgorithms.ZipFile
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    mktempdir() do dir
        path = joinpath(dir, "ala_fit.3mf")
        export_ses_3mf(sys, path; density=2.0, max_size_mm=100)
        # Read vertex positions back out of the 3MF and check max extent.
        xmin = ymin = zmin = Inf
        xmax = ymax = zmax = -Inf
        r = ZipFile.Reader(path)
        for f in r.files
            endswith(f.name, "3dmodel.model") || continue
            xml = String(read(f))
            for m in eachmatch(r"<vertex x=\"([-\d.e+]+)\" y=\"([-\d.e+]+)\" z=\"([-\d.e+]+)\" />", xml)
                x = parse(Float64, m.captures[1])
                y = parse(Float64, m.captures[2])
                z = parse(Float64, m.captures[3])
                xmin = min(xmin, x); xmax = max(xmax, x)
                ymin = min(ymin, y); ymax = max(ymax, y)
                zmin = min(zmin, z); zmax = max(zmax, z)
            end
        end
        close(r)
        # The mesh's longest extent is slightly larger than the atom-bbox
        # extent (the probe inflates the surface ~1.5 Å × scale on each
        # side). `max_size_mm` is anchored to the atom bbox, not the mesh
        # bbox, so we expect ~100 mm + 2 · 1.5 · (100/atom_extent) on the
        # longest axis. For AlaAla (atom extent ≈ 8.7 Å, scale ≈ 11.5 mm/Å)
        # the upper bound is ≈ 134 mm; allow [95, 145] for headroom across
        # Float32/Float64 and across Julia versions.
        max_ext = max(xmax - xmin, ymax - ymin, zmax - zmin)
        @test 95 <= max_ext <= 145
    end
end

@testitem "export_ses_3mf by_chain emits one 3MF object per chain" begin
    ZipFile = BiochemicalAlgorithms.ZipFile
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    nchains = length(chains(sys))
    @test nchains >= 1
    mktempdir() do dir
        path = joinpath(dir, "ala_bychain.3mf")
        export_ses_3mf(sys, path; density=2.0, by_chain=true)
        @test isfile(path)
        # Count <object> entries in the 3MF model XML.
        n_objects = 0
        r = ZipFile.Reader(path)
        for f in r.files
            endswith(f.name, "3dmodel.model") || continue
            xml = String(read(f))
            n_objects = length(collect(eachmatch(r"<object ", xml)))
            break
        end
        close(r)
        @test n_objects == nchains
    end
end

@testitem "construction_kit parts sit on z=0 and carry face_colors" begin
    sys = System{Float64}()
    mol = Molecule(sys; name="water")
    ch = Chain(mol); f = Fragment(ch, 1; name="WAT")
    oa = Atom(f, 1, Elements.O; r=Vector3{Float64}(0, 0, 0))
    ha = Atom(f, 2, Elements.H; r=Vector3{Float64}(0.96, 0, 0))
    hb = Atom(f, 3, Elements.H; r=Vector3{Float64}(-0.24, 0.93, 0))
    Bond(sys, oa.idx, ha.idx, BondOrder.Single)
    Bond(sys, oa.idx, hb.idx, BondOrder.Single)
    assign_radii!(sys)
    parts = construction_kit(sys; scale=20, flat_base_mm=0)
    # Parts must not overlap on the build plate (grid layout).
    function xy_bbox(p)
        xs = [Float64(v[1]) for v in p.mesh.position]
        ys = [Float64(v[2]) for v in p.mesh.position]
        (minimum(xs), maximum(xs), minimum(ys), maximum(ys))
    end
    boxes = [xy_bbox(p) for p in parts]
    for p in parts
        zs = [v[3] for v in p.mesh.position]
        # z-min == 0 (sits on build plate, since flat_base_mm=0).
        @test isapprox(minimum(zs), 0; atol=1e-8)
        # Every part has uniform face_colors set.
        @test p.face_colors !== nothing
        @test length(p.face_colors) == length(p.mesh.faces)
        @test length(unique(p.face_colors)) == 1
    end

    # With the default flat_base_mm=0.5, every part should sit `0.5 mm`
    # below the build plate so the slicer clips a small flat disk at z=0.
    parts2 = construction_kit(sys; scale=20)
    for p in parts2
        zmin = minimum(v[3] for v in p.mesh.position)
        @test isapprox(zmin, -0.5; atol=1e-8)
    end
    # No pair of parts has overlapping xy bounding boxes.
    for i in eachindex(boxes), j in (i+1):length(boxes)
        (x0i, x1i, y0i, y1i) = boxes[i]
        (x0j, x1j, y0j, y1j) = boxes[j]
        overlap = (x0i < x1j) && (x0j < x1i) && (y0i < y1j) && (y0j < y1i)
        @test !overlap
    end
end

@testitem "construction_kit parts are watertight (χ=2, no boundary, no non-manifold)" begin
    function audit(mesh)
        V = length(mesh.position); F = length(mesh.faces)
        ec = Dict{Tuple{Int,Int}, Int}()
        for f in mesh.faces
            for (a, b) in ((Int(f[1]), Int(f[2])), (Int(f[2]), Int(f[3])), (Int(f[3]), Int(f[1])))
                k = a < b ? (a, b) : (b, a)
                ec[k] = get(ec, k, 0) + 1
            end
        end
        E = length(ec)
        boundary = count(v -> v == 1, values(ec))
        non_manifold = count(v -> v > 2, values(ec))
        (V = V, E = E, F = F, χ = V - E + F, boundary = boundary, non_manifold = non_manifold)
    end

    # Methane: 4 sockets on C + 1 socket on each H + 4 bonds with pegs.
    sys = System{Float64}()
    mol = Molecule(sys; name="methane")
    ch = Chain(mol); f = Fragment(ch, 1; name="MET")
    ca = Atom(f, 1, Elements.C; r=Vector3{Float64}(0, 0, 0))
    positions = [Vector3{Float64}(1,1,1)*0.629, Vector3{Float64}(-1,-1,1)*0.629,
                 Vector3{Float64}(1,-1,-1)*0.629, Vector3{Float64}(-1,1,-1)*0.629]
    for (i, p) in enumerate(positions)
        h = Atom(f, i+1, Elements.H; r=p)
        Bond(sys, ca.idx, h.idx, BondOrder.Single)
    end
    assign_radii!(sys)
    parts = construction_kit(sys; scale=20)
    @test length(parts) == 9
    for p in parts
        s = audit(p.mesh)
        # χ = 2 ⇒ closed orientable genus-0 surface (= topological sphere)
        @test s.χ == 2
        @test s.boundary == 0
        @test s.non_manifold == 0
    end

    # Sanity: the C atom (4 sockets) should have more vertices than the H
    # atoms (1 socket each) — proving sockets are actually being cut.
    c_part = parts[findfirst(p -> endswith(p.name, "-C"), parts)]
    h_part = parts[findfirst(p -> endswith(p.name, "-H"), parts)]
    @test length(c_part.mesh.position) > length(h_part.mesh.position)
end

@testitem "export_ses_3mf reorient=true reduces z extent" begin
    ZipFile = BiochemicalAlgorithms.ZipFile
    sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), Float64)
    assign_radii!(sys)
    ses = compute_ses(sys; probe_radius=1.5)

    function zextent(path)
        zmin = Inf; zmax = -Inf
        r = ZipFile.Reader(path)
        for f in r.files
            endswith(f.name, "3dmodel.model") || continue
            xml = String(read(f))
            for m in eachmatch(r"<vertex [^/]*z=\"([-\d.e+]+)\" />", xml)
                z = parse(Float64, m.captures[1])
                zmin = min(zmin, z); zmax = max(zmax, z)
            end
        end
        close(r)
        zmax - zmin
    end

    mktempdir() do dir
        unori = joinpath(dir, "unori.3mf")
        ori   = joinpath(dir, "ori.3mf")
        export_ses_3mf(ses, sys, unori; density=2.0, reorient=false)
        export_ses_3mf(ses, sys, ori;   density=2.0, reorient=true)
        # PCA-aligned mesh: longest axis on X, shortest on Z. AlaAla is
        # elongated along the peptide backbone, so reorient should shave at
        # least 20% off the print height.
        z_unori = zextent(unori)
        z_ori   = zextent(ori)
        @test z_ori < z_unori
        @test z_ori <= z_unori * 0.85
    end
end
