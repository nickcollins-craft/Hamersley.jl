#=
Interpolations.

Fields declare a Ferrite interpolation family and an order. This file holds the single
gate deciding which families are supported (for now only Lagrange), and the resolution of
a family against a reference shape into a concrete Ferrite interpolation. Supporting a
further family (Serendipity, discontinuous Lagrange, Raviart–Thomas, …) means adding a
method to `is_supported_family`, plus a test, and nothing else in the package should need
to change, because:

- conformity (``H^1``, ``L^2``, ``H(\mathrm{div})``, ``H(\mathrm{curl})``) is read from
  Ferrite rather than declared;
- only scalar families are vectorised, so natively vector-valued families pass through;
- assembly goes through Ferrite's cell values and never assumes that degrees of freedom
  are nodal values.
=#


# ---------------------------------------------------------------------------------------- #
# --------------------------------------- The gate --------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    SUPPORTED_FAMILIES

The Ferrite interpolation families that Hamersley currently supports.
"""
const SUPPORTED_FAMILIES = (Ferrite.Lagrange,)


"""
    is_supported_family(spec)

Return `true` when the interpolation family of `spec` is supported. This is the only
place that decides it.
"""
is_supported_family(::FamilySpec) = false
is_supported_family(::FamilySpec{Ferrite.Lagrange}) = true


"""
    supported_family_names()

Return the names of the supported families, for error messages.
"""
supported_family_names() = join(string.(nameof.(SUPPORTED_FAMILIES)), ", ")


"""
    unsupported_family_message(field)

Return the message reporting that the family of `field` is not supported, listing the
families that are. The formulation verification and [`check_supported_family`](@ref) share it, so that the
two always say the same thing.
"""
function unsupported_family_message(field::Field)
    family = interpolation_family(field_interpolation(field))
    message =
        "field $(field_name(field)) uses $family, which is not yet supported, " *
            "please choose from: $(supported_family_names())"
    return message
end


"""
    check_supported_family(field)

Throw an `ArgumentError` naming `field`, its family and the supported families if the
family of `field` is not supported.
"""
function check_supported_family(field::Field)
    # Unlike the formulation verification, which collects issues, we throw at once: nothing that follows can
    # proceed without a supported family.
    if !is_supported_family(field_interpolation(field))
        throw(ArgumentError(unsupported_family_message(field)))
    end
    return nothing
end


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- Conformity -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    REPRESENTATIVE_REFERENCE_SHAPES

Reference shapes tried, in order, when we need an instance of a family only to ask Ferrite
for its conformity (which depends on the family, not the shape).
"""
const REPRESENTATIVE_REFERENCE_SHAPES = (
    Ferrite.RefTriangle,
    Ferrite.RefQuadrilateral,
    Ferrite.RefTetrahedron,
    Ferrite.RefHexahedron,
    Ferrite.RefLine,
)


"""
    family_conformity(spec)

Return Ferrite's conformity (for example `Ferrite.H1Conformity()`) of the family of
`spec`.
"""
function family_conformity(spec::FamilySpec)
    family = interpolation_family(spec)
    order = interpolation_order(spec)

    # Conformity is a property of the family, so any reference shape the family can be
    # instantiated on will do.
    for reference_shape in REPRESENTATIVE_REFERENCE_SHAPES
        interpolation = try
            family{reference_shape,order}()
        catch
            nothing
        end
        if interpolation !== nothing
            return Ferrite.conformity(interpolation)
        end
    end
    throw(
        ArgumentError("cannot instantiate $family of order $order on any reference shape"),
    )
end


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- Resolution -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    resolve_interpolation(field, reference_shape, dimension)

Return the concrete Ferrite interpolation of `field` on cells of `reference_shape` in a
mesh of spatial dimension `dimension`.

The family is resolved against the reference shape here, per subdomain, because Ferrite
interpolations are parametrised by it (and some families do not exist on every shape).
Scalar families are vectorised for vector fields; natively vector-valued families are
used as they are.
"""
function resolve_interpolation(
    field::Field,
    reference_shape::Type{<:Ferrite.AbstractRefShape},
    dimension::Integer,
)
    # The gate comes first, so that unsupported families fail with a clear message.
    check_supported_family(field)
    spec = field_interpolation(field)
    family = interpolation_family(spec)
    base_interpolation = family{reference_shape,interpolation_order(spec)}()

    # Not every family and order exists on every reference shape; Ferrite only finds out
    # when asked for the basis, so we ask now.
    try
        Ferrite.getnbasefunctions(base_interpolation)
    catch error_value
        if error_value isa MethodError
            throw(
                ArgumentError(
                    "$family of order $(interpolation_order(spec)) is not available on " *
                        "$reference_shape (field $(field_name(field)))",
                ),
            )
        else
            rethrow()
        end
    end

    interpolation = vectorise(field_shape(field), base_interpolation, dimension)
    return interpolation
end


"""
    vectorise(shape, interpolation, dimension)

Adapt `interpolation` to a field of shape `shape`: scalar interpolations are vectorised
for vector fields (``\\mathrm{ip}^d``), vector interpolations are kept as they are.
"""
vectorise(::ScalarField, interpolation::Ferrite.ScalarInterpolation, ::Integer) =
    interpolation
vectorise(::VectorField, interpolation::Ferrite.ScalarInterpolation, dimension::Integer) =
    interpolation^dimension
vectorise(::VectorField, interpolation::Ferrite.VectorInterpolation, ::Integer) =
    interpolation

function vectorise(::ScalarField, interpolation::Ferrite.VectorInterpolation, ::Integer)
    throw(
        ArgumentError("a scalar field cannot use the vector interpolation $interpolation"),
    )
end

function vectorise(
    shape::Union{SymmetricTensorField,TensorField},
    interpolation::Ferrite.Interpolation,
    dimension::Integer,
)
    not_yet_implemented("interpolation of $(nameof(typeof(shape))) fields")
    return nothing
end
