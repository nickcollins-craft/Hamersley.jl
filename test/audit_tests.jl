@testset "Problem audit" begin
    u = Field(:u, VectorField(), H1Space(); configuration=ReferenceConfiguration())
    T = Field(:T, ScalarField(), H1Space(); configuration=ReferenceConfiguration())

    mechanics = StressResponse(
        stress=FirstPiolaStress(),
        conjugate_to=DeformationGradient(:u),
        measure=ReferenceVolume(),
        dependencies=(:u, :T),
        law=PlaceholderLaw(:thermoplasticity),
    )
    heat = FluxResponse(
        flux=SpatialHeatFlux(:T),
        driven_by=SpatialGradient(:T),
        measure=CurrentVolume(),
        dependencies=(:T, :u),
        law=PlaceholderLaw(:fourier_finite_deformation),
    )

    problem = Problem(
        fields=(u, T),
        equations=(
            Equation(:momentum, :u, TestGradientTerm(mechanics)),
            Equation(:heat, :T, TestGradientTerm(heat)),
        ),
        responses=(mechanics, heat),
        boundary_conditions=(
            NaturalBC(:u, ValueTrace(); on=:right, quantity=ReferenceTraction()),
            EssentialBC(:T, ValueTrace(); on=:left, value=PlaceholderValue(:Tbar)),
        ),
    )
    @test passed(audit(problem))
    @test occursin("fields: u, T", explain(problem))

    undeclared_dependency = StressResponse(
        stress=FirstPiolaStress(),
        conjugate_to=DeformationGradient(:u),
        measure=ReferenceVolume(),
        dependencies=(:u, :p),
        law=PlaceholderLaw(:bad_dependency),
    )
    bad_problem = Problem(fields=(u,), responses=(undeclared_dependency,))
    @test !passed(audit(bad_problem))

    bad_equation = Problem(fields=(u,), equations=(Equation(:mass, :p),))
    @test !passed(audit(bad_equation))
end
