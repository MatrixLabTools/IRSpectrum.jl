using IRSpectrum
using Test


@testset "Read CP2K traj file" begin
    d_old  = read_cp2k_dipoles(joinpath(@__DIR__, "testdata", "dipoles-old.traj"))
    d_new  = read_cp2k_dipoles(joinpath(@__DIR__, "testdata", "dipoles-new.traj"))
    @test size(d_old) == (3, 5)
    @test size(d_new) == (3,4)

end