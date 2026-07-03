abstract type AbstractFieldShape end
struct ScalarField <: AbstractFieldShape end
struct VectorField <: AbstractFieldShape end
struct TensorField <: AbstractFieldShape end

abstract type AbstractFunctionSpace end
struct H1Space <: AbstractFunctionSpace end
struct H2Space <: AbstractFunctionSpace end
struct L2Space <: AbstractFunctionSpace end

abstract type AbstractFieldRole end
struct UnknownField <: AbstractFieldRole end
struct TestField <: AbstractFieldRole end
struct InternalVariableField <: AbstractFieldRole end

struct Field{Name,S<:AbstractFieldShape,F<:AbstractFunctionSpace,C<:AbstractConfiguration,R<:AbstractFieldRole}
    shape::S
    space::F
    configuration::C
    role::R
end

function Field(
    name::Symbol,
    shape::AbstractFieldShape,
    space::AbstractFunctionSpace;
    configuration::AbstractConfiguration=ReferenceConfiguration(),
    role::AbstractFieldRole=UnknownField(),
)
    return Field{name,typeof(shape),typeof(space),typeof(configuration),typeof(role)}(
        shape,
        space,
        configuration,
        role,
    )
end

field_name(::Field{Name}) where {Name} = Name
field_shape(field::Field) = field.shape
field_space(field::Field) = field.space
field_configuration(field::Field) = field.configuration
field_role(field::Field) = field.role

supports_strong_trace(::H1Space, ::Any) = false
supports_strong_trace(::H1Space, ::ValueTrace) = true
supports_strong_trace(::H2Space, ::Any) = true
supports_strong_trace(::L2Space, ::Any) = false
