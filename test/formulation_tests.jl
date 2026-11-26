#=
Tests of the formulation verification: the cases that motivated the restructure (small-strain
elasticity, finite-strain PK1 and PK2, and a pulled-back thermo-mechanical problem), and
the mistakes it must catch.
=#


"""
    has_issue(report, fragment)

Return `true` when some issue message in `report` contains `fragment`.
"""
has_issue(report, fragment) = any(contains(fragment), Hamersley.issue_messages(report))


@testset "Formulation verification" begin
    u = Field(:u, VectorField(), Lagrange, 1)
    F = DeformationGradient(:u)
    finite = FiniteStrain(motion=:u)

    @testset "Valid formulations" begin
        # Small-strain linear elasticity: rejected by the previous structure.
        @test passed(verification(small_strain_elasticity()))

        # Finite-strain neo-Hookean, ∫ P : δF dv_0.
        neo_hookean = Constitutive(
            NeoHookean(shear_modulus=1.0, first_lame_parameter=2.0),
            (F=F,),
        )
        pk1 = Problem(;
            fields=(u,),
            equations=(
                Equation(
                    :momentum,
                    :u,
                    Term(
                        Variation(F),
                        output(neo_hookean, :P),
                        ReferenceVolume(),
                    ),
                ),
            ),
            kinematics=finite,
        )
        @test passed(verification(pk1))

        # Finite strains with S and δE, ∫ S : δE dv_0.
        E = GreenLagrangeStrain(:u)
        saint_venant = Constitutive(
            PlaceholderLaw(
                :saint_venant_kirchhoff;
                inputs=(E=GreenLagrangeStrainTensor(),),
                outputs=(S=SecondPiolaStress(),),
            ),
            (E=E,),
        )
        pk2 = Problem(;
            fields=(u,),
            equations=(
                Equation(
                    :momentum,
                    :u,
                    Term(
                        Variation(E),
                        output(saint_venant, :S),
                        ReferenceVolume(),
                    ),
                ),
            ),
            kinematics=finite,
        )
        @test passed(verification(pk2))

        # Finite-strain thermo-mechanics with a pulled-back spatial heat flux.
        @test passed(verification(finite_strain_thermomechanics()))
    end

    @testset "Cauchy stress over the reference volume" begin
        cauchy_law = Constitutive(
            PlaceholderLaw(
                :spatial;
                inputs=(F=DeformationGradientTensor(),),
                outputs=(σ=CauchyStress(),),
            ),
            (F=F,),
        )
        term = Term(Variation(F), output(cauchy_law, :σ), ReferenceVolume())
        problem = Problem(;
            fields=(u,),
            equations=(Equation(:momentum, :u, term),),
            kinematics=finite,
        )
        report = verification(problem)
        @test !passed(report)
        @test has_issue(report, "not work-conjugate")
        @test has_issue(report, "wrap the flux in PullBack")
    end

    @testset "Infinitesimal strain under finite strains" begin
        problem = Problem(;
            fields=(u,),
            equations=Hamersley.problem_equations(small_strain_elasticity()),
            kinematics=finite,
        )
        @test has_issue(verification(problem), "only meaningful under SmallStrain")
    end

    @testset "Test sides" begin
        # A term whose variation acts on a different field than the equation's.
        T = Field(:T, ScalarField(), Lagrange, 1; represents=Temperature())
        source = Constitutive(
            PlaceholderLaw(
                :source;
                inputs=(T=Temperature(),),
                outputs=(r=ReferenceVolumeDensity(ScalarField()),),
            ),
            (T=Value(:T),),
        )
        wrong_field = Problem(;
            fields=(u, T),
            equations=(
                Equation(
                    :heat,
                    :T,
                    Term(
                        Variation(Value(:u)),
                        output(source, :r),
                        ReferenceVolume(),
                    ),
                ),
            ),
            kinematics=SmallStrain(),
        )
        @test has_issue(verification(wrong_field), "is tested by T")

        # Undeclared fields and test fields.
        undeclared = Problem(;
            fields=(u,),
            equations=(Equation(:heat, :T),),
            kinematics=SmallStrain(),
        )
        @test has_issue(verification(undeclared), "undeclared test field T")
    end

    @testset "Motion field" begin
        problem = Problem(; fields=(u,), kinematics=FiniteStrain(motion=:v))
        @test has_issue(verification(problem), "motion field v is not declared")
    end

    @testset "Unknown outputs" begin
        neo_hookean = Constitutive(
            NeoHookean(shear_modulus=1.0, first_lame_parameter=2.0),
            (F=F,),
        )
        problem = Problem(;
            fields=(u,),
            equations=(
                Equation(
                    :momentum,
                    :u,
                    Term(
                        Variation(F),
                        output(neo_hookean, :S),
                        ReferenceVolume(),
                    ),
                ),
            ),
            kinematics=finite,
        )
        @test has_issue(verification(problem), "please choose from: P")
    end

    @testset "Throwing on error" begin
        problem = Problem(; fields=(u,), kinematics=FiniteStrain(motion=:v))
        @test_throws ArgumentError verification(problem; throw_on_error=true)
    end
end
