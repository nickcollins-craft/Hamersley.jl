#=
Tests of the core vocabulary: kinematic assumptions, configuration collapse, the
pull-back tables and the conjugacy table.
=#


@testset "Core vocabulary" begin
    small = SmallStrain()
    finite = FiniteStrain(motion = :u)

    @testset "Kinematic assumptions" begin
        # A problem must declare its assumption; it is never guessed.
        @test_throws ArgumentError Problem(fields = ())
        @test Hamersley.motion_field(finite) === :u
        @test Hamersley.motion_field(small) === nothing
    end

    @testset "Configuration collapse" begin
        # Under small strains the current configuration is identified with the reference.
        @test Hamersley.same_configuration(small, CauchyStress(), ReferenceVolume())
        @test !Hamersley.same_configuration(finite, CauchyStress(), ReferenceVolume())
        @test Hamersley.same_configuration(finite, CauchyStress(), CurrentVolume())

        # The Kirchhoff stress is spatial but measured per unit reference volume.
        @test Hamersley.same_configuration(finite, KirchhoffStress(), ReferenceVolume())

        # Configuration-independent quantities are compatible with anything.
        @test Hamersley.same_configuration(finite, Temperature(), CurrentVolume())
    end

    @testset "Pull-backs and push-forwards" begin
        @test Hamersley.pulled_back(CauchyStress()) === FirstPiolaStress()
        @test Hamersley.pulled_back(SpatialFlux()) === MaterialFlux()
        @test Hamersley.pushed_forward(MaterialGradient(ScalarField())) ===
              SpatialGradient(ScalarField())
        @test Hamersley.pulled_back(SecondPiolaStress()) === nothing
        @test Hamersley.transformation_rule(SpatialFlux()) isa Hamersley.PiolaVectorRule
        @test Hamersley.transformation_rule(SpatialTraction()) isa Hamersley.NansonRule
    end

    @testset "Work conjugacy" begin
        are_conjugate = Hamersley.are_conjugate

        # Small strains: ε : σ, and every stress measure collapses onto σ.
        @test are_conjugate(small, InfinitesimalStrainTensor(), CauchyStress())
        @test are_conjugate(small, DeformationGradientTensor(), CauchyStress())
        @test are_conjugate(small, GreenLagrangeStrainTensor(), FirstPiolaStress())

        # Finite strains: δF : P and δE : S, but not δF : σ.
        @test are_conjugate(finite, DeformationGradientTensor(), FirstPiolaStress())
        @test are_conjugate(finite, GreenLagrangeStrainTensor(), SecondPiolaStress())
        @test !are_conjugate(finite, DeformationGradientTensor(), CauchyStress())

        # The table is symmetric, which is what stress-based mixed formulations need.
        @test are_conjugate(small, CauchyStress(), InfinitesimalStrainTensor())

        # A rate has the conjugates of the quantity itself.
        storage = Hamersley.RateOf(ReferenceVolumeDensity(ScalarField()))
        @test are_conjugate(finite, Temperature(), storage)
    end
end
