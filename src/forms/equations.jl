struct Equation{Terms<:Tuple}
    name::Symbol
    test_field::Symbol
    terms::Terms
end

function Equation(name::Symbol, test_field::Symbol, terms::AbstractWeakTerm...)
    return Equation{typeof(terms)}(name, test_field, terms)
end

equation_name(equation::Equation) = equation.name
test_field(equation::Equation) = equation.test_field
equation_terms(equation::Equation) = equation.terms
