abstract type AbstractOperator end

struct ValueOperator{FieldName} <: AbstractOperator end
struct ReferenceGradient{FieldName} <: AbstractOperator end
struct SpatialGradient{FieldName} <: AbstractOperator end
struct SymmetricGradient{FieldName,C<:AbstractConfiguration} <: AbstractOperator
    configuration::C
end
struct DeformationGradient{FieldName} <: AbstractOperator end
struct Divergence{FieldName,C<:AbstractConfiguration} <: AbstractOperator
    configuration::C
end

ValueOperator(field::Symbol) = ValueOperator{field}()
ReferenceGradient(field::Symbol) = ReferenceGradient{field}()
SpatialGradient(field::Symbol) = SpatialGradient{field}()
SymmetricGradient(field::Symbol; configuration::AbstractConfiguration=ReferenceConfiguration()) =
    SymmetricGradient{field,typeof(configuration)}(configuration)
DeformationGradient(field::Symbol) = DeformationGradient{field}()
Divergence(field::Symbol; configuration::AbstractConfiguration=ReferenceConfiguration()) =
    Divergence{field,typeof(configuration)}(configuration)

operator_field(::ValueOperator{Name}) where {Name} = Name
operator_field(::ReferenceGradient{Name}) where {Name} = Name
operator_field(::SpatialGradient{Name}) where {Name} = Name
operator_field(::SymmetricGradient{Name}) where {Name} = Name
operator_field(::DeformationGradient{Name}) where {Name} = Name
operator_field(::Divergence{Name}) where {Name} = Name

operator_configuration(::ValueOperator) = nothing
operator_configuration(::ReferenceGradient) = ReferenceConfiguration()
operator_configuration(::SpatialGradient) = CurrentConfiguration()
operator_configuration(op::SymmetricGradient) = op.configuration
operator_configuration(::DeformationGradient) = ReferenceConfiguration()
operator_configuration(op::Divergence) = op.configuration
