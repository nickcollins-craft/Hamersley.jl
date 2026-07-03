@testset "Core metadata" begin
    u = Field(:u, VectorField(), H1Space(); configuration=ReferenceConfiguration())
    T = Field(:T, ScalarField(), H1Space(); configuration=ReferenceConfiguration())

    @test Hamersley.field_name(u) === :u
    @test Hamersley.field_name(T) === :T
    @test Hamersley.operator_field(DeformationGradient(:u)) === :u
    @test Hamersley.operator_configuration(ReferenceGradient(:T)) isa ReferenceConfiguration
    @test Hamersley.operator_configuration(SpatialGradient(:T)) isa CurrentConfiguration
    @test Hamersley.requires_higher_continuity(NormalGradientTrace())
end
