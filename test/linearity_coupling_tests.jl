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

    @testset "Follower loads" begin
        follower = finite_strain_thermomechanics(traction=SpatialTraction())
        bc = first(follower.boundary_conditions)
        @test Hamersley.is_follower_load(bc, follower.kinematics)
        @test !Hamersley.is_follower_load(bc, SmallStrain())
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
