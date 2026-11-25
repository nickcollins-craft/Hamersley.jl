#=
Tests of the rate factors that time integrators supply to the linearisation.
=#


@testset "Time integrators" begin
    Δt = 0.1

    @test Hamersley.rate_factor(QuasiStatic(), 1, Δt) == 0.0
    @test Hamersley.rate_factor(BackwardEuler(), 1, Δt) ≈ 1/Δt
    @test_throws ArgumentError Hamersley.rate_factor(BackwardEuler(), 2, Δt)

    # Average-acceleration Newmark: γ/(βΔt) and 1/(βΔt²).
    newmark = Newmark(beta = 1/4, gamma = 1/2)
    @test Hamersley.rate_factor(newmark, 1, Δt) ≈ (1/2)/((1/4)*Δt)
    @test Hamersley.rate_factor(newmark, 2, Δt) ≈ 1/((1/4)*Δt^2)
end
