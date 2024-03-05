@testitem "Ring Perception" begin
    system = load_sdfile(ball_data_path("../test/data/rings_test.sdf"))

    num_ring_atoms = [6, 6, 6, 6, 6, 13, 0, 6, 4]

    detected_num_ring_atoms = map(sum, map(is_ring_atom, molecules(system)))
    @test num_ring_atoms == detected_num_ring_atoms

    num_sssrs      = [1, 1, 1, 1, 1, 3, 0, 1, 3]
    num_sssr_atoms = [[6], [6], [6], [6], [6], [6, 6, 3], [], [6], [3, 3, 3]]

    detected_sssrs = map(find_sssr, molecules(system))
    @test map(length, detected_sssrs) == num_sssrs
    @test [map(length, sssr) for sssr in detected_sssrs] == num_sssr_atoms
end