#=
Tests of the interpolation gate, conformity, resolution against reference shapes, and the
grid checks.
=#


@testset "Interpolations" begin
    @testset "The gate" begin
        lagrange = Field(:u, VectorField(), Lagrange, 2)
        serendipity = Field(:u, VectorField(), Ferrite.Serendipity, 2)
        @test Hamersley.is_supported_family(lagrange.interpolation)
        @test !Hamersley.is_supported_family(serendipity.interpolation)
        @test_throws ArgumentError Hamersley.check_supported_family(serendipity)

        # The formulation verification reports unsupported families too, listing what is supported.
        report = verification(Problem(; fields=(serendipity,), kinematics=SmallStrain()))
        messages = Hamersley.issue_messages(report)
        @test any(contains("please choose from: Lagrange"), messages)
    end

    @testset "Conformity is read from Ferrite" begin
        conformity(family, order) =
            Hamersley.family_conformity(Hamersley.FamilySpec(family, order))
        @test conformity(Lagrange, 1) isa Ferrite.H1Conformity
        @test conformity(Ferrite.RaviartThomas, 1) isa Ferrite.HdivConformity
        @test conformity(Ferrite.DiscontinuousLagrange, 0) isa Ferrite.L2Conformity
    end

    @testset "Resolution" begin
        u = Field(:u, VectorField(), Lagrange, 2)
        T = Field(:T, ScalarField(), Lagrange, 1)

        # Scalar families are vectorised for vector fields only.
        @test Hamersley.resolve_interpolation(u, Ferrite.RefTriangle, 2) ==
            Lagrange{Ferrite.RefTriangle,2}()^2
        @test Hamersley.resolve_interpolation(T, Ferrite.RefHexahedron, 3) ==
            Lagrange{Ferrite.RefHexahedron,1}()
    end

    @testset "Grids" begin
        grid = Ferrite.generate_grid(Ferrite.Triangle, (2, 2))
        problem = small_strain_elasticity()
        @test isempty(Hamersley.grid_issues(grid, problem))

        # Facet sets the grid lacks are reported.
        base = small_strain_elasticity()
        missing_set = Problem(;
            fields=base.fields,
            equations=base.equations,
            boundary_conditions=(
                NaturalBC(:u; on=:lateral, quantity=ReferenceTraction()),
            ),
            kinematics=SmallStrain(),
        )
        @test_throws ArgumentError Hamersley.validate_grid(grid, missing_set)

        # Compilation runs both checks before reaching what is not yet written.
        @test_throws NotYetImplemented compile(problem, grid)
    end
end
