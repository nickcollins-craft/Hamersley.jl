#=
The nodes of a formulation's dataflow graph.

A Hamersley formulation is a small directed graph:

    fields → operators → tagged quantities → laws → tagged quantities → weak terms

This file defines the abstract node type, the operators that act directly on a field
(value, gradient, divergence, second gradient), and the two wrappers that every weak form
needs: `Variation`, which marks the test side of a term, and `Rate`, a time derivative
whose relation to the unknowns is supplied by the time integrator.

Every node answers four questions, each by a generic function:

- `quantity(node, fields)`: which quantity tag its value carries;
- `node_fields(node, assumption)`: which fields its value depends on;
- `is_affine(node, assumption)`: whether its value is affine in those fields;
- `evaluate` and `derivative`: its value and linearisation at a point (not yet written).

Linearisation belongs to every node, not just to constitutive laws: under finite strains
the geometric and load stiffnesses come from kinematic and transformation nodes.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Abstract nodes ------------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractNode

Supertype of every node of a formulation's dataflow graph.
"""
abstract type AbstractNode end


"""
    AbstractFieldOperator

Nodes that apply a differential operator directly to a single field.
"""
abstract type AbstractFieldOperator <: AbstractNode end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Field operators ----------------------------------- #
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
field is not declared (the audit reports undeclared fields separately).
"""
function field_quantity(tag_of_field, operator::AbstractNode, fields::Tuple)
    field = lookup_field(fields, operand_field(operator))
    if field === nothing
        return nothing
    else
        return tag_of_field(field)
    end
end


node_fields(operator::AbstractFieldOperator, ::AbstractKinematicAssumption) =
    [operand_field(operator)]

# The field operators are linear in their field.
is_affine(::AbstractFieldOperator, ::AbstractKinematicAssumption) = true

subnodes(::AbstractFieldOperator) = ()


# ---------------------------------------------------------------------------------------- #
# ------------------------------------- Variations --------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    Variation(node)

The variation ``δ(\\cdot)`` of `node` in the direction of the test function: the test side
of a weak term. `Variation(DeformationGradient(:u))` is ``δF = ∇_0 δu``, and
`Variation(GreenLagrangeStrain(:u))` is ``δE = \\operatorname{sym}(F^T ∇_0 δu)``, which
depends on the solution and therefore contributes to the tangent.

The node being varied must depend on exactly one field, which is the test field of the
equation the term belongs to.

# Fields
- `node::N`: the node being varied.
"""
struct Variation{N<:AbstractNode} <: AbstractNode
    node::N
end


quantity(variation::Variation, fields::Tuple) = quantity(variation.node, fields)

# The variation is constant (independent of the unknowns) exactly when the node is affine,
# so it contributes to the tangent through the node's own field dependencies.
node_fields(variation::Variation, assumption::AbstractKinematicAssumption) =
    node_fields(variation.node, assumption)
is_affine(variation::Variation, assumption::AbstractKinematicAssumption) =
    is_affine(variation.node, assumption)

subnodes(variation::Variation) = (variation.node,)


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------- Rates ----------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    Rate(node; order = 1)

The time derivative of order `order` of `node`. A law receives a rate as an ordinary
input and differentiates with respect to it; the time integrator supplies the factor
``∂\\dot{x}/∂x_{n+1}`` that relates it to the current unknowns (for backward Euler,
``1/Δt``).

A rate may wrap a field operator (the rate of the field, then the operator) or a derived
quantity (the rate of the quantity itself, as in ``\\dot{J} ≈ (J_{n+1} - J_n)/Δt``). The
two agree as ``Δt → 0`` but not at finite ``Δt``, so the choice is part of the
formulation.

# Fields
- `node::N`: the node whose rate is taken.
- `order::Int`: the order of the time derivative (1 for ``\\dot{x}``, 2 for ``\\ddot{x}``).
"""
struct Rate{N<:AbstractNode} <: AbstractNode
    node::N
    order::Int
end


function Rate(node::AbstractNode; order::Integer = 1)
    # Only first and second time derivatives appear in the formulations we support.
    if !(order == 1 || order == 2)
        throw(ArgumentError("Rate order must be 1 or 2, got $order"))
    end
    rate = Rate(node, Int(order))
    return rate
end


function quantity(rate::Rate, fields::Tuple)
    # A rate of an unknown quantity is itself unknown.
    inner = quantity(rate.node, fields)
    if inner === nothing
        return nothing
    else
        return RateOf(inner)
    end
end

node_fields(rate::Rate, assumption::AbstractKinematicAssumption) =
    node_fields(rate.node, assumption)
is_affine(rate::Rate, assumption::AbstractKinematicAssumption) =
    is_affine(rate.node, assumption)

subnodes(rate::Rate) = (rate.node,)


# ---------------------------------------------------------------------------------------- #
# ------------------------------- Evaluation and linearisation --------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    evaluate(node, point_state)

Evaluate `node` at a quadrature point, given the interpolated field values and gradients
in `point_state`.
"""
function evaluate(node::AbstractNode, point_state)
    not_yet_implemented("evaluation of $(nameof(typeof(node))) nodes")
    return nothing
end


"""
    derivative(node, field, point_state)

Return the derivative of `node` with respect to the field called `field` at a quadrature
point. The compile step chains these derivatives from each weak term back to the fields,
which is how every block of a multiphysics tangent is formed.
"""
function derivative(node::AbstractNode, field::Symbol, point_state)
    not_yet_implemented("linearisation of $(nameof(typeof(node))) nodes")
    return nothing
end
