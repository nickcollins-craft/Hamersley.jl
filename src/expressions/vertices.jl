#=
The vertices of a formulation's dataflow graph.

A Hamersley formulation is a small directed graph:

    fields → operators → tagged quantities → constitutive laws → tagged quantities → weak terms

This file defines the abstract vertex (vertex is to be understood in the graph-theoretic
sense) type, the operators that act directly on a field (e.g. value, gradient, divergence,
second gradient and so on), and the two wrappers that every weak form needs: `Variation`,
which marks the test side of a term, and `Rate`, a time derivative whose relation to the
unknowns is supplied by the time integrator.

Every vertex carries four pieces of information, each conveyed by using a generic function:

- `quantity(vertex, fields)`: which quantity tag its value carries;
- `vertex_fields(vertex, assumption)`: which fields its value depends on;
- `is_affine(vertex, assumption)`: whether its value is affine in those fields;
- `evaluate` and `derivative`: its value and linearisation at a point.

Because under finite strains the derivation needs to propagate backwards through every
vertex in the graph, linearisation is relevant to every vertex, not just to constitutive
laws: under finite strains the geometric and load stiffnesses come from the kinematic and
transformation vertices, respectively. For the trivial vertices, is_affine tells the solver
that the tangent is constant and only needs to be factorised once
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Abstract vertices ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractVertex

Supertype of every vertex of a formulation's dataflow graph.
"""
abstract type AbstractVertex end


"""
    AbstractFieldOperator

Vertices that apply a differential operator directly to a single field.
"""
abstract type AbstractFieldOperator <: AbstractVertex end


"""
    quantity(vertex, fields)

Return the quantity tag carried by the value of `vertex`, given the tuple of declared
`fields`, or `nothing` when it cannot be known (for example because `vertex` refers to an
undeclared field, which the formulation verification reports separately).
"""
function quantity end


"""
    vertex_fields(vertex, assumption)

Return a vector of the names of the fields that the value of `vertex` depends on under
`assumption`. This is what the coupling matrix and the sparsity pattern are built from.
"""
function vertex_fields end


"""
    subvertices(vertex)

Return a tuple of the vertices directly beneath `vertex` in the graph.
"""
function subvertices end


"""
    operand_field(vertex)

Return the name of the field that the field operator or kinematic vertex `vertex` acts on.
"""
function operand_field end


# ---------------------------------------------------------------------------------------- #
# ----------------------------------- Field operators ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    Value(field)

The point value of `field`. Its quantity tag is whatever the field declares it represents,
which is how interpolated stresses and pressures enter mixed formulations.

# Fields
- `field::Symbol`: the name of the field.
"""
struct Value <: AbstractFieldOperator
    field::Symbol
end


"""
    ReferenceGradient(field)

The gradient ``∇_0`` of `field` with respect to the material coordinates.

# Fields
- `field::Symbol`: the name of the field.
"""
struct ReferenceGradient <: AbstractFieldOperator
    field::Symbol
end


"""
    Divergence(field)

The material divergence ``\\operatorname{Div}`` of the vector field `field`.

# Fields
- `field::Symbol`: the name of the field.
"""
struct Divergence <: AbstractFieldOperator
    field::Symbol
end


"""
    SecondGradient(field)

The second gradient ``∇_0 ∇_0`` of `field`. This is part of the structure for
second-gradient continua; with ``C^0`` (Lagrange) elements it cannot be discretised
directly, and such models are written as micromorphic formulations instead.

# Fields
- `field::Symbol`: the name of the field.
"""
struct SecondGradient <: AbstractFieldOperator
    field::Symbol
end


operand_field(operator::AbstractFieldOperator) = operator.field

quantity(operator::Value, fields::Tuple) =
    field_quantity(field_value_quantity, operator, fields)
quantity(operator::ReferenceGradient, fields::Tuple) =
    field_quantity(field -> MaterialGradient(field_shape(field)), operator, fields)
quantity(operator::Divergence, fields::Tuple) =
    field_quantity(_ -> MaterialDivergence(), operator, fields)
quantity(operator::SecondGradient, fields::Tuple) =
    field_quantity(_ -> MaterialSecondGradient(), operator, fields)


"""
    field_quantity(tag_of_field, operator, fields)

Apply `tag_of_field` to the field that `operator` acts on, or return `nothing` if that
field is not declared (the formulation verification reports undeclared fields separately).
"""
function field_quantity(tag_of_field, operator::AbstractVertex, fields::Tuple)
    # An undeclared field has no shape or meaning, so it gives no tag.
    field = lookup_field(fields, operand_field(operator))
    if field === nothing
        return nothing
    else
        return tag_of_field(field)
    end
end


"""
    wrapped_quantity(tag_of_inner, wrapper, fields)

Apply `tag_of_inner` to the quantity tag of the vertex that `wrapper` wraps (a rate, a
pull-back or a push-forward), or return `nothing` if that inner tag is unknown.
"""
function wrapped_quantity(tag_of_inner, wrapper::AbstractVertex, fields::Tuple)
    # A rate or transform of an unknown quantity is itself unknown.
    inner = quantity(wrapper.vertex, fields)
    if inner === nothing
        return nothing
    else
        return tag_of_inner(inner)
    end
end


vertex_fields(operator::AbstractFieldOperator, ::AbstractKinematicAssumption) =
    [operand_field(operator)]

# The field operators are linear in their field.
is_affine(::AbstractFieldOperator, ::AbstractKinematicAssumption) = true

subvertices(::AbstractFieldOperator) = ()


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- Variations -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    Variation(vertex)

The variation ``δ(\\cdot)`` of `vertex` in the direction of the test function: the test side
of a weak term. `Variation(DeformationGradient(:u))` is ``δF = ∇_0 δu``, and
`Variation(GreenLagrangeStrain(:u))` is ``δE = \\operatorname{sym}(F^T ∇_0 δu)``, which
depends on the solution and therefore contributes to the tangent.

The vertex being varied must depend on exactly one field, which is the test field of the
equation the term belongs to.

# Fields
- `vertex::N`: the vertex being varied.
"""
struct Variation{N<:AbstractVertex} <: AbstractVertex
    vertex::N
end


quantity(variation::Variation, fields::Tuple) = quantity(variation.vertex, fields)

# The variation is constant (independent of the unknowns) exactly when the vertex is affine,
# so it contributes to the tangent through the vertex's own field dependencies.
vertex_fields(variation::Variation, assumption::AbstractKinematicAssumption) =
    vertex_fields(variation.vertex, assumption)
is_affine(variation::Variation, assumption::AbstractKinematicAssumption) =
    is_affine(variation.vertex, assumption)

subvertices(variation::Variation) = (variation.vertex,)


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------- Rates ----------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    Rate(vertex; order = 1)

The time derivative of order `order` of `vertex`. A law receives a rate as an ordinary
input and differentiates with respect to it; the time integrator supplies the factor
``∂\\dot{x}/∂x_{n+1}`` that relates it to the current unknowns (for backward Euler,
``1/Δt``).

A rate may wrap a field operator (the rate of the field, then the operator) or a derived
quantity (the rate of the quantity itself, as in ``\\dot{J} ≈ (J_{n+1} - J_n)/Δt``). The
two agree as ``Δt → 0`` but not at finite ``Δt``, so the choice is part of the
formulation.

# Fields
- `vertex::N`: the vertex whose rate is taken.
- `order::Int`: the order of the time derivative (1 for ``\\dot{x}``, 2 for ``\\ddot{x}``).
"""
struct Rate{N<:AbstractVertex} <: AbstractVertex
    vertex::N
    order::Int
end


function Rate(vertex::AbstractVertex; order::Integer=1)
    # Only first and second time derivatives appear in the formulations we support.
    if !(order == 1 || order == 2)
        throw(ArgumentError("Rate order must be 1 or 2, got $order"))
    end
    rate = Rate(vertex, Int(order))
    return rate
end


# A rate carries the tag of the quantity it is the rate of, wrapped in `RateOf`.
quantity(rate::Rate, fields::Tuple) = wrapped_quantity(RateOf, rate, fields)

vertex_fields(rate::Rate, assumption::AbstractKinematicAssumption) =
    vertex_fields(rate.vertex, assumption)
is_affine(rate::Rate, assumption::AbstractKinematicAssumption) =
    is_affine(rate.vertex, assumption)

subvertices(rate::Rate) = (rate.vertex,)


# ---------------------------------------------------------------------------------------- #
# ----------------------------- Evaluation and linearisation ----------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    evaluate(vertex, point_state)

Evaluate `vertex` at a quadrature point, given the interpolated field values and gradients
in `point_state`.
"""
function evaluate(vertex::AbstractVertex, point_state)
    not_yet_implemented("evaluation of $(nameof(typeof(vertex))) vertices")
    return nothing
end


"""
    derivative(vertex, field, point_state)

Return the derivative of `vertex` with respect to the field called `field` at a quadrature
point. The compile step chains these derivatives from each weak term back to the fields,
which is how every block of a multiphysics tangent is formed.
"""
function derivative(vertex::AbstractVertex, field::Symbol, point_state)
    not_yet_implemented("linearisation of $(nameof(typeof(vertex))) vertices")
    return nothing
end
