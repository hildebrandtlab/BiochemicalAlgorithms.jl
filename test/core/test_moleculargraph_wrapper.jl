@testitem "Convert to SDFMolGraph" begin
    using MolecularGraph

    # Simple system
    for T in [Float32, Float64]
        sys = load_pdb(ball_data_path("../test/data/AlaAla.pdb"), T)

        fdb = FragmentDB{T}()
        normalize_names!(sys, fdb)
        reconstruct_fragments!(sys, fdb)
        build_bonds!(sys, fdb)

        g = convert(SDFMolGraph, sys)
        @test g isa SDFMolGraph
        aidx_dict = g.gprops[:atom_idx]

        @test length(g.vprops) == natoms(sys)
        for (i, at) in enumerate(atoms(sys))
            @test g.vprops[i].symbol == Symbol(at.element)
            @test g.vprops[i].coords ≈ at.r
            @test aidx_dict[i] == at.idx
        end

        @test length(g.eprops) == nbonds(sys)
        bds_sys = Set((b.a1, b.a2, Int(b.order)) for b in bonds(sys))
        bds_g   = Set((aidx_dict[e.first.src], aidx_dict[e.first.dst], e.second.order) for e in g.eprops)
        @test length(bds_sys ∩ bds_g) == nbonds(sys)
    end

    # Non-unique atom numbers
    # https://github.com/hildebrandtlab/BiochemicalAlgorithms.jl/issues/142
    for T in [Float32, Float64]
        sys = System{T}()
        Atom(sys, 1, Elements.H)
        Atom(sys, 1, Elements.C)

        g = convert(SDFMolGraph, sys)
        @test g isa SDFMolGraph
        aidx_dict = g.gprops[:atom_idx]

        @test length(g.vprops) == natoms(sys)
        for (i, at) in enumerate(atoms(sys))
            @test g.vprops[i].symbol == Symbol(at.element)
            @test g.vprops[i].coords ≈ at.r
            @test aidx_dict[i] == at.idx
        end
    end
end
