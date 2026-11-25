#=
Tests of the derivation of boundary traces from the weak form, of natural conditions that
infer their trace, of strict boundary checking, and of the report given by `explain`.
=#


@testset "Boundary traces" begin
    @testset "Derivation from test operators" begin
        problem = small_strain_elasticity()
        @test Hamersley.derived_traces(problem, :u) == [ValueTrace()]

        # A second-gradient term leaves the value, the normal gradient and edge values.
        u = Field(:u, VectorField(), Lagrange, 2)
        hyperstress = Constitutive(
            PlaceholderLaw(
                :second_gradient;
                inputs = (G = MaterialSecondGradient(),),
                outputs = (T = HyperStress(),),
            ),
            (G = SecondGradient(:u),),
        )
        second_gradient = Problem(;
            fields = (u,),
            equations = (
                Equation(
                    :momentum,
                    :u,
                    Term(
                        Variation(SecondGradient(:u)),
                        output(hyperstress, :T),
                        ReferenceVolume(),
                    ),
                ),
            ),
            boundary_conditions = (NaturalBC(:u; on = :top, quantity = DoubleForce()),),
            kinematics = SmallStrain(),
        )
        @test Hamersley.derived_traces(second_gradient, :u) ==
              [ValueTrace(), NormalGradientTrace(), EdgeValueTrace()]

        # The double force infers its trace, which is no longer ambiguous.
        @test passed(audit(second_gradient))

        # Strict mode insists on every non-standard trace being dealt with.
        strict = Problem(;
            fields = second_gradient.fields,
            equations = second_gradient.equations,
            boundary_conditions = (
                EssentialBC(:u; on = :bottom, value = PlaceholderValue(:u_bar)),
            ),
            kinematics = SmallStrain(),
            strict_boundaries = true,
        )
        report = audit(strict)
        messages = Hamersley.issue_messages(report)
        @test any(contains("NormalGradientTrace() of field u on bottom"), messages)
    end

    @testset "Traces that do not exist" begin
        # Prescribing a double force on a first-gradient field is rejected.
        base = small_strain_elasticity()
        problem = Problem(;
            fields = base.fields,
            equations = base.equations,
            boundary_conditions = (NaturalBC(:u; on = :top, quantity = DoubleForce()),),
            kinematics = SmallStrain(),
        )
        messages = Hamersley.issue_messages(audit(problem))
        @test any(contains("does not appear in the weak form"), messages)

        # Prescribing the normal gradient of a Lagrange field cannot be done strongly.
        strong_gradient = EssentialBC(
            :u;
            on = :top,
            value = PlaceholderValue(:g),
            trace = NormalGradientTrace(),
        )
        u = Field(:u, VectorField(), Lagrange, 1)
        @test Hamersley.imposition(strong_gradient, u) isa Hamersley.WeakImposition
    end

    @testset "Explain" begin
        text = explain(small_strain_elasticity())
        @test occursin("Problem under SmallStrain()", text)
        @test occursin("left: u ValueTrace() — Prescribed", text)
        @test occursin("right: u ValueTrace() — Loaded", text)
    end
end
