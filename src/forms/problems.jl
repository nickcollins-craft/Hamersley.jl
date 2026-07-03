struct FormulationMetadata
    supported_traces::Dict{Symbol,Tuple}
    higher_order_fields::Tuple{Vararg{Symbol}}
end

function FormulationMetadata(;
    supported_traces::Dict{Symbol,<:Tuple}=Dict{Symbol,Tuple}(),
    higher_order_fields::Tuple{Vararg{Symbol}}=(),
)
    return FormulationMetadata(Dict{Symbol,Tuple}(supported_traces), higher_order_fields)
end

struct Problem{Fields<:Tuple,Equations<:Tuple,Responses<:Tuple,BCs<:Tuple}
    fields::Fields
    equations::Equations
    responses::Responses
    boundary_conditions::BCs
    metadata::FormulationMetadata
end

function Problem(;
    fields::Tuple,
    equations::Tuple=(),
    responses::Tuple=(),
    boundary_conditions::Tuple=(),
    metadata::FormulationMetadata=FormulationMetadata(),
)
    return Problem{typeof(fields),typeof(equations),typeof(responses),typeof(boundary_conditions)}(
        fields,
        equations,
        responses,
        boundary_conditions,
        metadata,
    )
end

problem_fields(problem::Problem) = problem.fields
problem_equations(problem::Problem) = problem.equations
problem_responses(problem::Problem) = problem.responses
problem_boundary_conditions(problem::Problem) = problem.boundary_conditions
problem_metadata(problem::Problem) = problem.metadata

function field_names(problem::Problem)
    names = Symbol[]
    for field in problem_fields(problem)
        push!(names, field_name(field))
    end
    return names
end

function equation_names(problem::Problem)
    names = Symbol[]
    for equation in problem_equations(problem)
        push!(names, equation_name(equation))
    end
    return names
end

function field_by_name(problem::Problem, name::Symbol)
    for field in problem_fields(problem)
        field_name(field) === name && return field
    end
    return nothing
end

function supported_trace_types(problem::Problem, field::Symbol)
    explicit = get(problem_metadata(problem).supported_traces, field, nothing)
    explicit !== nothing && return explicit
    return (ValueTrace,)
end

function supports_trace(problem::Problem, field::Symbol, trace::AbstractTrace)
    return typeof(trace) in supported_trace_types(problem, field)
end

function is_higher_order_field(problem::Problem, field::Symbol)
    return field in problem_metadata(problem).higher_order_fields
end
