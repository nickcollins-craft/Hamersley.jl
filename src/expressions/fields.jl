#=
Fields.

A field is an unknown of the problem, interpolated by finite elements. It declares its
shape, its interpolation (a Ferrite family and an order, resolved against the mesh's
reference shapes later), and optionally the physical quantity it represents. That last
declaration is what makes mixed formulations ordinary: an interpolated pressure or stress
field carries the same quantity tag as a constitutive output, and the formulation verification treats the two
alike.

Fields do not carry a configuration. Under the total Lagrangian approach every field is
parametrised over the reference configuration; configurations belong to quantities and
measures.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Interpolations ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    FamilySpec(family, order)

A finite element family and polynomial order, not yet tied to a reference shape.

The family is a Ferrite interpolation family such as `Ferrite.Lagrange`, stored as a type
parameter so that support for it can be decided by dispatch. It is resolved to a concrete
Ferrite interpolation (for example `Lagrange{RefTriangle, 2}()`) once the mesh is known;
see `discretisation/interpolations.jl`.

# Fields
- `order::Int`: the polynomial order of the interpolation.
"""
struct FamilySpec{F}
    order::Int
end


function FamilySpec(family::Type{<:Ferrite.Interpolation}, order::Integer)
    # We guard against meaningless orders here rather than at resolution time.
    if order < 0
        throw(ArgumentError("interpolation order must be non-negative, got $order"))
    end
    spec = FamilySpec{family}(Int(order))
    return spec
end


"""
    interpolation_family(spec)

Return the Ferrite family (for example `Ferrite.Lagrange`) of `spec`.
"""
interpolation_family(::FamilySpec{F}) where {F} = F


"""
    interpolation_order(spec)

Return the polynomial order of `spec`.
"""
interpolation_order(spec::FamilySpec) = spec.order


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------- Fields ---------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    Field(name, shape, family, order; represents = nothing, element_local = false)

An unknown field of the problem.

# Fields
- `name::Symbol`: the symbol of the field (for example `:u` or `:T`).
- `shape::S`: its [`AbstractFieldShape`](@ref).
- `interpolation::I`: its [`FamilySpec`](@ref).
- `represents::Q`: the quantity tag the field's value carries (for example
  `Pressure()` or `CauchyStress()`), or `nothing` for a generic
  [`FieldValue`](@ref).
- `element_local::Bool`: whether the field's degrees of freedom belong to single
  elements, so that they can be statically condensed. Condensation is not yet
  implemented, and the formulation verification reports such fields.

# Examples
```julia
u = Field(:u, VectorField(), Lagrange, 2)
T = Field(:T, ScalarField(), Lagrange, 1; represents = Temperature())
```
"""
struct Field{S<:AbstractFieldShape,I<:FamilySpec,Q}
    name::Symbol
    shape::S
    interpolation::I
    represents::Q
    element_local::Bool
end


function Field(
    name::Symbol,
    shape::AbstractFieldShape,
    family::Type{<:Ferrite.Interpolation},
    order::Integer;
    represents::Union{AbstractQuantity,Nothing}=nothing,
    element_local::Bool=false,
)
    spec = FamilySpec(family, order)
    field = Field(name, shape, spec, represents, element_local)
    return field
end


field_name(field::Field) = field.name
field_shape(field::Field) = field.shape
field_interpolation(field::Field) = field.interpolation


"""
    field_value_quantity(field)

Return the quantity tag carried by the point value of `field`: the tag it declares it
represents, or a generic [`FieldValue`](@ref) of its shape.
"""
function field_value_quantity(field::Field)
    # A field that declares no physical meaning is a generic value of its shape.
    if field.represents === nothing
        return FieldValue(field.shape)
    else
        return field.represents
    end
end


"""
    lookup_field(fields, name)

Return the field called `name` in the tuple `fields`, or `nothing` if there is none.
"""
function lookup_field(fields::Tuple, name::Symbol)
    # Problems have a handful of fields, so a linear search is the clearest choice.
    for field in fields
        if field_name(field) === name
            return field
        end
    end
    return nothing
end
