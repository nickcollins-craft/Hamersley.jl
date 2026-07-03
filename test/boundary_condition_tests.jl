@testset "Boundary condition traces" begin
    u = Field(:u, VectorField(), H1Space(); configuration=ReferenceConfiguration())
    metadata = FormulationMetadata(
        supported_traces=Dict(:u => (ValueTrace, NormalGradientTrace)),
        higher_order_fields=(:u,),
    )

    second_gradient_problem = Problem(
        fields=(u,),
        boundary_conditions=(
            NaturalBC(:u, NormalGradientTrace(); on=:top, quantity=DoubleForce()),
        ),
        metadata=metadata,
    )
    @test passed(audit(second_gradient_problem))

    ambiguous = Problem(
        fields=(u,),
        boundary_conditions=(
            NaturalBC(:u; on=:top, quantity=DoubleForce()),
        ),
        metadata=metadata,
    )
    @test !passed(audit(ambiguous))

    unsupported_trace = Problem(
        fields=(u,),
        boundary_conditions=(
            NaturalBC(:u, TangentialComponentTrace(); on=:top, quantity=DoubleForce()),
        ),
        metadata=metadata,
    )
    @test !passed(audit(unsupported_trace))

    strong_normal_gradient = Problem(
        fields=(u,),
        boundary_conditions=(
            EssentialBC(:u, NormalGradientTrace(); on=:top, value=PlaceholderValue(:dubar_dn)),
        ),
        metadata=metadata,
    )
    @test !passed(audit(strong_normal_gradient))
end
