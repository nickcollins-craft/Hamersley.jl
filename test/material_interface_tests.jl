#=
Tests of the constitutive law interface and the library declarations.
=#


@testset "Material interface" begin
    @testset "Library declarations" begin
        elasticity =
            IsotropicLinearElasticity(first_lame_parameter=1, shear_modulus=2.0)
        @test elasticity.λ === 1.0
        @test Hamersley.law_inputs(elasticity) == (ε=InfinitesimalStrainTensor(),)
        @test Hamersley.law_outputs(elasticity) == (σ=CauchyStress(),)
        @test Hamersley.is_affine(elasticity)
        @test Hamersley.tangent_is_symmetric(elasticity)

        neo_hookean = NeoHookean(shear_modulus=1.0, first_lame_parameter=2.0)
        @test Hamersley.law_inputs(neo_hookean) == (F=DeformationGradientTensor(),)
        @test !Hamersley.is_affine(neo_hookean)
        @test Hamersley.tangent_is_symmetric(neo_hookean)

        # By default every output depends on every input.
        @test Hamersley.output_dependencies(neo_hookean) == (P=(:F,),)
    end

    @testset "Placeholder laws" begin
        law = PlaceholderLaw(
            :coupled;
            inputs=(F=DeformationGradientTensor(), T=Temperature()),
            outputs=(P=FirstPiolaStress(),),
        )
        @test Hamersley.output_dependencies(law) == (P=(:F, :T),)
        @test Hamersley.law_name(law) === :coupled
        @test_throws NotYetImplemented Hamersley.constitutive_update(law, (;), nothing, 1.0)
    end

    @testset "Bindings" begin
        u = Field(:u, VectorField(), Lagrange, 1)
        neo_hookean = NeoHookean(shear_modulus=1.0, first_lame_parameter=2.0)

        # A binding with the right tag has no issues.
        good = Constitutive(neo_hookean, (F=DeformationGradient(:u),))
        @test isempty(Hamersley.binding_issues(good, (u,)))

        # A binding with the wrong tag, or an unknown slot, is reported.
        wrong_tag = Constitutive(neo_hookean, (F=GreenLagrangeStrain(:u),))
        @test length(Hamersley.binding_issues(wrong_tag, (u,))) == 1
        unknown_slot =
            Constitutive(neo_hookean, (F=DeformationGradient(:u), E=Value(:u)))
        @test length(Hamersley.binding_issues(unknown_slot, (u,))) == 1
    end

    @testset "Local solves" begin
        @test_throws NotYetImplemented Hamersley.solve_local(identity, 0.0, nothing)
    end
end
