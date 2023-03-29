@testitem "Ring Perception" begin
    system = load_sdfile(ball_data_path("../test/data/rings_test.sdf"))

    num_ring_atoms = [6, 6, 6, 6, 6, 13, 0, 6, 4]

    detected_num_ring_atoms = map(length, map(is_ring_atom, molecules(system)))

    @test num_ring_atoms == detected_num_ring_atoms
end