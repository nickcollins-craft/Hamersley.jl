#=
Tests of the inference of linearity, the check of a user's linearity claim, the default
solver strategy, and the field coupling matrix.
=#


@testset "Linearity and coupling" begin
    @testset "Linearity" begin
        @test Hamersley.is_affine(small_strain_elasticity())
        @test passed(verification(small_strain_elasticity(linear=true)))
        @test Hamersley.default_strategy(small_strain_elasticity()) isa
            CachedFactorisationStrategy

        # Finite strains are not linear, so a claim that they are fails the formulation verification.
        thermomechanics = finite_strain_thermomechanics()
        @test !Hamersley.is_affine(thermomechanics)
        @test Hamersley.default_strategy(thermomechanics) isa NewtonStrategy
    end

    @testset "Opting out of linearity" begin
        # An affine formulation declared nonlinear is solved by Newton's method.
        opted_out = small_strain_elasticity(linear=false)
        @test Hamersley.is_affine(opted_out)
        @test passed(verification(opted_out))
        @test Hamersley.default_strategy(opted_out) isa NewtonStrategy
    end

    @testset "Viscoelasticity" begin
        # Linear Maxwell viscoelasticity, a law of ε and its rate carrying σ_n as state.
        u = Field(:u, VectorField(), Lagrange, 1)
        ε = InfinitesimalStrain(:u)
        function maxwell(; affine, linear=nothing)
            law = Constitutive(
                PlaceholderLaw(
                    :maxwell;
                    inputs=(
                        ε=InfinitesimalStrainTensor(),
                        ε_dot=RateOf(InfinitesimalStrainTensor()),
                    ),
                    outputs=(σ=CauchyStress(),),
                    affine=affine,
                    initial_state=(σ_n=zeros(Float64, 6),),
                ),
                (ε=ε, ε_dot=Rate(ε)),
            )
            problem = Problem(;
                fields=(u,),
                equations=(
                    Equation(
                        :momentum,
                        :u,
                        Term(Variation(ε), output(law, :σ), ReferenceVolume()),
                    ),
                ),
                kinematics=SmallStrain(),
                linear=linear,
            )
            return problem
        end

        # A law declared affine makes the problem affine, and its state is still seen.
        affine_maxwell = maxwell(affine=true)
        @test passed(verification(affine_maxwell))
        @test Hamersley.is_affine(affine_maxwell)
        @test Hamersley.has_internal_state(affine_maxwell)
        @test Hamersley.default_strategy(affine_maxwell) isa CachedFactorisationStrategy
        @test Hamersley.default_strategy(maxwell(affine=true, linear=false)) isa
            NewtonStrategy

        # Without the declaration, laws are assumed nonlinear.
        @test !Hamersley.is_affine(maxwell(affine=false))
        @test !Hamersley.has_internal_state(small_strain_elasticity())
    end

    @testset "Follower loads" begin
        # A follower traction becomes a term over the current boundary, which depends on
        # the motion under finite strains and not under small strains.
        follower = finite_strain_thermomechanics(traction=SpatialTraction())
        bc = first(Hamersley.phase_conditions(follower))
        term = Hamersley.condition_term(follower, bc)
        @test Hamersley.term_measure(term) isa CurrentBoundary
        @test !Hamersley.is_affine(term, follower.kinematics)
        @test Hamersley.is_affine(term, SmallStrain())
        @test Hamersley.term_dependencies(term, follower.kinematics) == [:u]

        # A dead traction becomes a term over the reference boundary.
        dead = finite_strain_thermomechanics()
        dead_term = Hamersley.condition_term(dead, first(Hamersley.phase_conditions(dead)))
        @test Hamersley.term_measure(dead_term) isa ReferenceBoundary
        @test Hamersley.is_affine(dead_term, dead.kinematics)
    end

    @testset "Current measures" begin
        # An affine momentum equation in u alone, and heat conduction (affine) plus a term
        # integrated over a given measure, so that only the measure can make the heat
        # equation nonlinear or couple it to u.
        u = Field(:u, VectorField(), Lagrange, 2)
        T = Field(:T, ScalarField(), Lagrange, 1; represents=Temperature())
        F = DeformationGradient(:u)
        linear_piola = Constitutive(
            PlaceholderLaw(
                :linear_piola;
                inputs=(F=DeformationGradientTensor(),),
                outputs=(P=FirstPiolaStress(),),
                affine=true,
            ),
            (F=F,),
        )
        momentum = Equation(
            :momentum,
            :u,
            Term(Variation(F), output(linear_piola, :P), ReferenceVolume()),
        )
        fourier = Constitutive(
            PlaceholderLaw(
                :fourier;
                inputs=(g=MaterialGradient(ScalarField()),),
                outputs=(q=MaterialFlux(),),
                affine=true,
            ),
            (g=ReferenceGradient(:T),),
        )
        function heat_with(term; kinematics=FiniteStrain(motion=:u), linear=nothing)
            conduction = Term(
                Variation(ReferenceGradient(:T)),
                output(fourier, :q),
                ReferenceVolume(),
            )
            problem = Problem(;
                fields=(u, T),
                equations=(momentum, Equation(:heat, :T, conduction, term)),
                kinematics=kinematics,
                linear=linear,
            )
            return problem
        end

        # A source r = c T per unit current volume: ∫ δT c T dv = ∫ δT c T J dv_0.
        source = Constitutive(
            PlaceholderLaw(
                :source;
                inputs=(T=Temperature(),),
                outputs=(r=SpatialVolumeDensity(ScalarField()),),
                affine=true,
            ),
            (T=Value(:T),),
        )
        current_source = Term(Variation(Value(:T)), output(source, :r), CurrentVolume())

        # The term is well formed, but under finite strains it depends on u through J.
        finite = heat_with(current_source)
        @test passed(verification(finite))
        @test !Hamersley.is_affine(finite)
        @test Hamersley.coupling_matrix(finite) == [true false; true true]
        @test has_issue(verification(heat_with(current_source; linear=true)), "not affine")

        # Under small strains J ≈ 1, so the same term is affine and uncoupled.
        small = heat_with(current_source; kinematics=SmallStrain())
        @test Hamersley.is_affine(small)
        @test Hamersley.coupling_matrix(small) == [true false; false true]

        # Over the reference volume the measure is fixed.
        reference_source = Constitutive(
            PlaceholderLaw(
                :source;
                inputs=(T=Temperature(),),
                outputs=(r=ReferenceVolumeDensity(ScalarField()),),
                affine=true,
            ),
            (T=Value(:T),),
        )
        reference = heat_with(
            Term(Variation(Value(:T)), output(reference_source, :r), ReferenceVolume()),
        )
        @test Hamersley.is_affine(reference)
        @test Hamersley.coupling_matrix(reference) == [true false; false true]

        # A convection term h T da over the current boundary picks up Nanson's formula.
        convection = heat_with(
            Term(Variation(Value(:T)), Value(:T), CurrentBoundary(); on=:right),
        )
        @test !Hamersley.is_affine(convection)
        @test Hamersley.coupling_matrix(convection) == [true false; true true]
    end

    @testset "Coupling" begin
        @test Hamersley.coupling_matrix(small_strain_elasticity()) == reshape([true], 1, 1)

        # Momentum depends on T through the law; heat depends on u through the transforms.
        thermomechanics = finite_strain_thermomechanics()
        @test Hamersley.coupling_matrix(thermomechanics) == [true true; true true]

        # A displacement–pressure saddle point: the pressure equation does not depend on p.
        u = Field(:u, VectorField(), Lagrange, 2)
        p = Field(:p, ScalarField(), Lagrange, 1; represents=Pressure())
        incompressibility = Equation(
            :incompressibility,
            :p,
            Term(Variation(Value(:p)), Divergence(:u), ReferenceVolume()),
        )
        saddle = Problem(;
            fields=(u, p),
            equations=(incompressibility,),
            kinematics=SmallStrain(),
        )
        @test passed(verification(saddle))
        @test Hamersley.coupling_matrix(saddle) == [false false; true false]
    end
end
