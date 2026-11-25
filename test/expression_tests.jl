#=
Tests of the dataflow graph: node tags, field dependencies, transforms, and affinity.
=#


@testset "Expressions" begin
    u = Field(:u, VectorField(), Lagrange, 1)
    p = Field(:p, ScalarField(), Lagrange, 1; represents = Pressure())
    fields = (u, p)
    small = SmallStrain()
    finite = FiniteStrain(motion = :u)

    @testset "Quantity tags" begin
        @test Hamersley.quantity(Value(:u), fields) === FieldValue(VectorField())
        @test Hamersley.quantity(Value(:p), fields) === Pressure()
        @test Hamersley.quantity(ReferenceGradient(:p), fields) ===
              MaterialGradient(ScalarField())
        @test Hamersley.quantity(DeformationGradient(:u), fields) ===
              DeformationGradientTensor()
        @test Hamersley.quantity(Rate(Value(:p)), fields) === Hamersley.RateOf(Pressure())

        # Undeclared fields give no tag (the audit reports them separately).
        @test Hamersley.quantity(Value(:missing), fields) === nothing
    end

    @testset "Transforms" begin
        pushed = PushForward(ReferenceGradient(:p))
        @test Hamersley.quantity(pushed, fields) === SpatialGradient(ScalarField())

        # Under finite strains a transform depends on the motion; under small strains not.
        @test Hamersley.node_fields(pushed, finite) == [:p, :u]
        @test Hamersley.node_fields(pushed, small) == [:p]
        @test !Hamersley.is_affine(pushed, finite)
        @test Hamersley.is_affine(pushed, small)
    end

    @testset "Kinematic affinity" begin
        @test Hamersley.is_affine(DeformationGradient(:u), finite)
        @test !Hamersley.is_affine(GreenLagrangeStrain(:u), finite)
        @test Hamersley.is_affine(GreenLagrangeStrain(:u), small)
        @test !Hamersley.is_affine(VolumeRatio(:u), finite)
    end

    @testset "Per-output dependencies" begin
        # The flux depends only on the gradient, so it never couples to the value.
        law = PlaceholderLaw(
            :transport;
            inputs = (g = MaterialGradient(ScalarField()), J = JacobianDeterminant()),
            outputs = (Q = MaterialFlux(), m = ReferenceVolumeDensity(ScalarField())),
            dependencies = (Q = (:g,), m = (:J,)),
        )
        binding = Constitutive(law, (g = ReferenceGradient(:p), J = VolumeRatio(:u)))
        @test Hamersley.node_fields(output(binding, :Q), finite) == [:p]
        @test Hamersley.node_fields(output(binding, :m), finite) == [:u]
    end

    @testset "Rates" begin
        @test_throws ArgumentError Rate(Value(:p); order = 3)
        @test Rate(Value(:p); order = 2).order == 2
    end

    @testset "Placeholders" begin
        @test_throws NotYetImplemented Hamersley.evaluate(Value(:u), nothing)
        @test_throws NotYetImplemented Hamersley.derivative(Value(:u), :u, nothing)
    end
end
