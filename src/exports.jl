#=
The public interface of Hamersley.

Only names a user types to build and run a problem are exported. Internals (traversal
helpers, formulation verification rules, transformation tables) stay reachable as `Hamersley.name`.
=#


# Errors.
export NotYetImplemented

# Configurations and measures.
export ReferenceConfiguration, CurrentConfiguration, IntermediateConfiguration
export ReferenceVolume, CurrentVolume, ReferenceBoundary, CurrentBoundary
export ReferenceInterface

# Kinematic assumptions.
export SmallStrain, FiniteStrain

# Field shapes.
export ScalarField, VectorField, SymmetricTensorField, TensorField

# Quantity tags.
export FieldValue, Temperature, Pressure
export MaterialGradient, SpatialGradient, MaterialDivergence, MaterialSecondGradient
export DeformationGradientTensor, GreenLagrangeStrainTensor, InfinitesimalStrainTensor
export JacobianDeterminant
export CauchyStress, KirchhoffStress, FirstPiolaStress, SecondPiolaStress, HyperStress
export MaterialFlux, SpatialFlux, ReferenceVolumeDensity, SpatialVolumeDensity
export ReferenceTraction, SpatialTraction, DoubleForce, EdgeForce
export ReferenceNormalFlux, SpatialNormalFlux, BoundaryPressure
export ReferencePressure, SpatialPressure
export RateOf

# Traces, regions and prescribed values.
export ValueTrace, NormalGradientTrace, NormalComponentTrace, TangentialComponentTrace
export EdgeValueTrace
export FacetSet, EdgeSet, VertexSet
export Constant, Prescribed, Tabulated

# Constitutive laws.
export PlaceholderLaw, IsotropicLinearElasticity, NeoHookean

# Fields and the graph. `Lagrange` is Ferrite's own family, re-exported for convenience.
export Field, Lagrange
export Value, ReferenceGradient, Divergence, SecondGradient, Variation, Rate
export DeformationGradient, GreenLagrangeStrain, InfinitesimalStrain, VolumeRatio
export PullBack, PushForward
export Constitutive, output

# Formulation.
export Term, Equation, EssentialBC, NaturalBC, RobinBC, Problem
export Phase, Schedule, AtTime, Never, Release
export verification, passed, explain

# Time integration and solution strategies.
export QuasiStatic, BackwardEuler, Newmark
export NewtonStrategy, CachedFactorisationStrategy, FullTangent, BlockDiagonalTangent

# Compilation and meshing.
export compile, Meshing, read_mesh
