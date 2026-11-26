#=
Problems.

A problem gathers the fields, the equations and the boundary conditions, together with
the one declaration that everything else depends on: the kinematic assumption. It also
records an optional user claim that the problem is linear, which the formulation verification checks once
against the structure of the formulation.
=#


"""
    Problem(; fields, equations = (), boundary_conditions = (), kinematics,
              linear = nothing, strict_boundaries = false, template = NoTemplate())

A formulation, ready to be verified and (eventually) compiled against a mesh.

# Fields
- `fields::Fs`: a tuple of [`Field`](@ref)s. Their order is the order of the fields in
  the degree-of-freedom handler and the coupling matrix.
- `equations::Es`: a tuple of [`Equation`](@ref)s, one per test field.
- `boundary_conditions::Bs`: a tuple of [`AbstractBoundaryCondition`](@ref)s.
- `kinematics::K`: the [`AbstractKinematicAssumption`](@ref), either `SmallStrain()` or
  `FiniteStrain(motion = :u)`.
- `linear::Union{Bool,Nothing}`: a user claim that the problem is (or is not) linear, or
  `nothing` to let the structure decide.
- `strict_boundaries::Bool`: whether every non-standard trace (such as the normal
  gradient of a second-gradient field) must be dealt with explicitly on every boundary.
- `template::T`: the [`AbstractFormulationTemplate`](@ref) that translates boundary
  conditions onto the discretised fields.
"""
struct Problem{
    Fs<:Tuple,
    Es<:Tuple,
    Bs<:Tuple,
    K<:AbstractKinematicAssumption,
    T<:AbstractFormulationTemplate,
}
    fields::Fs
    equations::Es
    boundary_conditions::Bs
    kinematics::K
    linear::Union{Bool,Nothing}
    strict_boundaries::Bool
    template::T
end


function Problem(;
    fields::Tuple,
    equations::Tuple=(),
    boundary_conditions::Tuple=(),
    kinematics::Union{AbstractKinematicAssumption,Nothing}=nothing,
    linear::Union{Bool,Nothing}=nothing,
    strict_boundaries::Bool=false,
    template::AbstractFormulationTemplate=NoTemplate(),
)
    # The kinematic assumption is never guessed: it decides which pairings are valid.
    if kinematics === nothing
        throw(
            ArgumentError(
                "a Problem must declare its kinematic assumption, please choose from: " *
                    "kinematics = SmallStrain(), kinematics = FiniteStrain(motion = :u)",
            ),
        )
    end

    problem = Problem(
        fields,
        equations,
        boundary_conditions,
        kinematics,
        linear,
        strict_boundaries,
        template,
    )
    return problem
end


problem_fields(problem::Problem) = problem.fields
problem_equations(problem::Problem) = problem.equations
problem_boundary_conditions(problem::Problem) = problem.boundary_conditions
kinematic_assumption(problem::Problem) = problem.kinematics


"""
    field_names(problem)

Return the names of the fields of `problem`, in declaration order.
"""
field_names(problem::Problem) = Symbol[field_name(field) for field in problem.fields]


"""
    equation_names(problem)

Return the names of the equations of `problem`, in declaration order.
"""
equation_names(problem::Problem) =
    Symbol[equation_name(equation) for equation in problem.equations]


"""
    field_by_name(problem, name)

Return the field of `problem` called `name`, or `nothing`.
"""
field_by_name(problem::Problem, name::Symbol) = lookup_field(problem.fields, name)


"""
    problem_terms(problem)

Return every weak term of `problem`, paired with the equation it belongs to, as a vector
of `(equation, term)` tuples.
"""
function problem_terms(problem::Problem)
    # We flatten the equations into one list, remembering which equation each term is in.
    pairs = Tuple{Equation,AbstractWeakTerm}[]
    for equation in problem.equations
        for term in equation_terms(equation)
            push!(pairs, (equation, term))
        end
    end
    return pairs
end


"""
    problem_constitutive(problem)

Return the distinct [`Constitutive`](@ref) bindings used anywhere in `problem`.
"""
function problem_constitutive(problem::Problem)
    bindings = Constitutive[]
    for (_, term) in problem_terms(problem)
        # A binding may feed several terms, so each one is only reported once.
        for vertex in (term_test(term), term_flux(term))
            for binding in collect_constitutive(vertex)
                push_unique!(bindings, binding)
            end
        end
    end
    return bindings
end
