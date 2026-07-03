abstract type AbstractQuantity end
abstract type AbstractStressMeasure <: AbstractQuantity end
abstract type AbstractStrainMeasure <: AbstractQuantity end
abstract type AbstractFluxQuantity <: AbstractQuantity end
abstract type AbstractSourceQuantity <: AbstractQuantity end
abstract type AbstractStorageQuantity <: AbstractQuantity end
abstract type AbstractBoundaryQuantity <: AbstractQuantity end

struct FirstPiolaStress <: AbstractStressMeasure end
struct SecondPiolaStress <: AbstractStressMeasure end
struct CauchyStress <: AbstractStressMeasure end
struct KirchhoffStress <: AbstractStressMeasure end

struct DeformationGradientQuantity <: AbstractStrainMeasure end
struct GreenLagrangeStrain <: AbstractStrainMeasure end
struct SpatialVelocityGradient <: AbstractStrainMeasure end

struct MaterialHeatFlux{FieldName} <: AbstractFluxQuantity end
struct SpatialHeatFlux{FieldName} <: AbstractFluxQuantity end
struct DarcyFlux{FieldName} <: AbstractFluxQuantity end

struct HeatSource <: AbstractSourceQuantity end
struct BodyForce <: AbstractSourceQuantity end
struct InternalEnergy <: AbstractStorageQuantity end
struct MassStorage <: AbstractStorageQuantity end

struct ReferenceTraction <: AbstractBoundaryQuantity end
struct SpatialTraction <: AbstractBoundaryQuantity end
struct DoubleForce <: AbstractBoundaryQuantity end
struct ReferenceFlux <: AbstractBoundaryQuantity end
struct SpatialFlux <: AbstractBoundaryQuantity end

struct PullBack{Q<:AbstractQuantity} <: AbstractQuantity
    quantity::Q
end

struct PlaceholderLaw
    name::Symbol
end

struct PlaceholderValue
    name::Symbol
end

MaterialHeatFlux(field::Symbol) = MaterialHeatFlux{field}()
SpatialHeatFlux(field::Symbol) = SpatialHeatFlux{field}()
DarcyFlux(field::Symbol) = DarcyFlux{field}()

quantity_field(::MaterialHeatFlux{Name}) where {Name} = Name
quantity_field(::SpatialHeatFlux{Name}) where {Name} = Name
quantity_field(::DarcyFlux{Name}) where {Name} = Name

quantity_configuration(::FirstPiolaStress) = ReferenceConfiguration()
quantity_configuration(::SecondPiolaStress) = ReferenceConfiguration()
quantity_configuration(::CauchyStress) = CurrentConfiguration()
quantity_configuration(::KirchhoffStress) = CurrentConfiguration()
quantity_configuration(::MaterialHeatFlux) = ReferenceConfiguration()
quantity_configuration(::SpatialHeatFlux) = CurrentConfiguration()
quantity_configuration(::DarcyFlux) = CurrentConfiguration()
quantity_configuration(::ReferenceTraction) = ReferenceConfiguration()
quantity_configuration(::SpatialTraction) = CurrentConfiguration()
quantity_configuration(::DoubleForce) = ReferenceConfiguration()
quantity_configuration(::ReferenceFlux) = ReferenceConfiguration()
quantity_configuration(::SpatialFlux) = CurrentConfiguration()
quantity_configuration(quantity::PullBack) = ReferenceConfiguration()

conjugate_configuration(::DeformationGradient) = ReferenceConfiguration()
conjugate_configuration(::ReferenceGradient) = ReferenceConfiguration()
conjugate_configuration(::SpatialGradient) = CurrentConfiguration()
conjugate_configuration(op::SymmetricGradient) = op.configuration
conjugate_configuration(::DeformationGradientQuantity) = ReferenceConfiguration()
conjugate_configuration(::GreenLagrangeStrain) = ReferenceConfiguration()
conjugate_configuration(::SpatialVelocityGradient) = CurrentConfiguration()

is_valid_stress_pair(::FirstPiolaStress, ::DeformationGradient) = true
is_valid_stress_pair(::SecondPiolaStress, ::GreenLagrangeStrain) = true
is_valid_stress_pair(::CauchyStress, ::SpatialVelocityGradient) = true
is_valid_stress_pair(::KirchhoffStress, ::SpatialVelocityGradient) = true
is_valid_stress_pair(::AbstractStressMeasure, ::Any) = false
