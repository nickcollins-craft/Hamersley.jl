#=
Equations.

An equation is the set of weak terms tested by one field. A multiphysics problem has one
equation per test field: balance of momentum tested by ``δu``, balance of energy tested by
``δθ``, and so on.
=#


"""
    Equation(name, test_field, terms...)

A weak equation: the sum of `terms`, tested by the field `test_field`, set to zero.

# Fields
- `name::Symbol`: a name for the equation (for example `:momentum`).
- `test_field::Symbol`: the field whose test function the terms are varied in.
- `terms::Terms`: a tuple of [`AbstractWeakTerm`](@ref)s.
"""
struct Equation{Terms<:Tuple}
    name::Symbol
    test_field::Symbol
    terms::Terms
end


function Equation(name::Symbol, test_field::Symbol, terms::AbstractWeakTerm...)
    equation = Equation(name, test_field, terms)
    return equation
end


equation_name(equation::Equation) = equation.name
test_field(equation::Equation) = equation.test_field
equation_terms(equation::Equation) = equation.terms
