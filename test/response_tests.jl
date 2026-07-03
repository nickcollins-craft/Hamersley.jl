@testset "Response metadata" begin
    mechanics = StressResponse(
        stress=FirstPiolaStress(),
        conjugate_to=DeformationGradient(:u),
        measure=ReferenceVolume(),
        dependencies=(:u, :T),
        law=PlaceholderLaw(:thermoplasticity),
    )
    @test isempty(Hamersley.response_issues(mechanics))
    @test Hamersley.response_dependencies(mechanics) == (:u, :T)

    invalid_stress = StressResponse(
        stress=CauchyStress(),
        conjugate_to=DeformationGradient(:u),
        measure=ReferenceVolume(),
        dependencies=(:u,),
        law=PlaceholderLaw(:bad_pair),
    )
    @test !isempty(Hamersley.response_issues(invalid_stress))

    finite_deformation_flux = FluxResponse(
        flux=SpatialHeatFlux(:T),
        driven_by=SpatialGradient(:T),
        measure=CurrentVolume(),
        dependencies=(:T, :u),
        law=PlaceholderLaw(:fourier_current),
    )
    @test isempty(Hamersley.response_issues(finite_deformation_flux))

    missing_pullback = FluxResponse(
        flux=SpatialHeatFlux(:T),
        driven_by=SpatialGradient(:T),
        measure=ReferenceVolume(),
        dependencies=(:T, :u),
        law=PlaceholderLaw(:missing_pullback),
    )
    @test !isempty(Hamersley.response_issues(missing_pullback))
end
