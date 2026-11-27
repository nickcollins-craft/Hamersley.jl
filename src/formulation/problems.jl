#=
Problems.

A problem gathers the fields, the equations and the boundary conditions (as a schedule of
phases, see `schedules.jl`), together with the one declaration that everything else
depends on: the kinematic assumption. It also records an optional user claim that the
problem is linear, which the formulation verification checks once against the structure of
the formulation.
=#


"""
    Problem(; fields, equations = (), boundary_conditions = (), schedule = nothing,
              kinematics, linear = nothing, strict_boundaries = false,
              template = NoTemplate())

A formulation, ready to be verified and (eventually) compiled against a mesh.

# Fields
- `fields::Fs`: a tuple of [`Field`](@ref)s. Their order is the order of the fields in
  the degree-of-freedom handler and the coupling matrix.
- `equations::Es`: a tuple of [`Equation`](@ref)s, one per test field.
- `schedule::S`: the [`Schedule`](@ref) of boundary conditions. Giving
  `boundary_conditions` instead (a tuple of [`AbstractBoundaryCondition`](@ref)s) makes a
  single phase that never ends; giving both is an error.
- `kinematics::K`: the [`AbstractKinematicAssumption`](@ref), either `SmallStrain()` or
  `FiniteStrain(motion = :u)`.
- `linear::Union{Bool,Nothing}`: a user claim that the problem is (or is not) linear, or
  `nothing` to let the structure decide. A claim of `true` is checked by the formulation
  verification; a claim of `false` is always honoured, and selects Newton's method even for
  an affine formulation.
- `strict_boundaries::Bool`: whether every non-standard trace (such as the normal
  gradient of a second-gradient field) must be dealt with explicitly on every boundary.
- `template::T`: the [`AbstractFormulationTemplate`](@ref) that translates boundary
  conditions onto the discretised fields.
"""
struct Problem{
    Fs<:Tuple,
    Es<:Tuple,
    S<:Schedule,
    K<:AbstractKinematicAssumption,
    T<:AbstractFormulationTemplate,
}
    fields::Fs
    equations::Es
    schedule::S
    kinematics::K
    linear::Union{Bool,Nothing}
    strict_boundaries::Bool
    template::T
end


function Problem(;
    fields::Tuple,
    equations::Tuple=(),
    boundary_conditions::Union{Tuple,Nothing}=nothing,
    schedule::Union{Schedule,Nothing}=nothing,
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

    # Boundary conditions are given either as one fixed set or as a schedule of phases.
    if boundary_conditions !== nothing && schedule !== nothing
        throw(
            ArgumentError(
                "a Problem takes boundary_conditions or schedule, not both; put the " *
                    "conditions in the first phase of the schedule",
            ),
        )
    end
    resolved_schedule = if schedule !== nothing
        schedule
    else
        single_phase_schedule(something(boundary_conditions, ()))
    end

    problem = Problem(
        fields,
        equations,
        resolved_schedule,
        kinematics,
        linear,
        strict_boundaries,
        template,
    )
    return problem
end


problem_fields(problem::Problem) = problem.fields
problem_equations(problem::Problem) = problem.equations
problem_schedule(problem::Problem) = problem.schedule
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


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------- Phases ---------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    problem_phases(problem)

Return the phases of the schedule of `problem`, in order.
"""
problem_phases(problem::Problem) = problem.schedule.phases


"""
    phase_index(problem, phase)

Return the position of the phase called `phase` in the schedule of `problem` (an integer
is returned as it is), throwing an `ArgumentError` that lists the phases otherwise.
"""
phase_index(::Problem, index::Integer) = Int(index)

function phase_index(problem::Problem, name::Symbol)
    names = Symbol[phase_name(phase) for phase in problem_phases(problem)]
    index = findfirst(==(name), names)
    if index === nothing
        valid = join(String[string(phase) for phase in names], ", ")
        throw(ArgumentError("the problem has no phase $name, please choose from: $valid"))
    end
    return index
end


"""
    phase_conditions(problem, phase = 1)

Return the boundary conditions that apply in `phase` (a name or a position) of `problem`,
after inheritance, replacement and releases.
"""
phase_conditions(problem::Problem, phase=1) =
    resolved_phases(problem.schedule)[phase_index(problem, phase)]


"""
    all_conditions(problem)

Return every boundary condition that applies in some phase of `problem`, once each.
"""
function all_conditions(problem::Problem)
    conditions = AbstractBoundaryCondition[]
    for conditions_of_phase in resolved_phases(problem.schedule)
        for bc in conditions_of_phase
            push_unique!(conditions, bc)
        end
    end
    return conditions
end


"""
    written_conditions(problem)

Return every boundary condition as written in some phase of `problem` (before inherited
conditions are restricted), once each. These are the conditions checked one by one.
"""
function written_conditions(problem::Problem)
    conditions = AbstractBoundaryCondition[]
    for phase in problem_phases(problem)
        for entry in phase.conditions
            if entry isa AbstractBoundaryCondition
                push_unique!(conditions, entry)
            end
        end
    end
    return conditions
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------- Weak terms --------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    weak_form_terms(problem)

Return the weak terms written in the equations of `problem`, paired with the equation they
belong to, as a vector of `(equation, term)` tuples.
"""
function weak_form_terms(problem::Problem)
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
    problem_terms(problem)
    problem_terms(problem, phase)

Return every weak term of `problem`, paired with the equation it belongs to: the terms
written in the equations, and the boundary terms that natural and Robin conditions become
(see [`condition_terms`](@ref)), in every phase or in the one phase `phase`.
"""
problem_terms(problem::Problem) =
    vcat(weak_form_terms(problem), condition_terms(problem, all_conditions(problem)))

problem_terms(problem::Problem, phase) = vcat(
    weak_form_terms(problem),
    condition_terms(problem, phase_conditions(problem, phase)),
)


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
