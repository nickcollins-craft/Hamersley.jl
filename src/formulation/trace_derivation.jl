#=
Derivation of boundary traces from the weak form.

Integrating a weak form by parts leaves one boundary term per trace of the test function,
and which traces appear depends only on the operators acting on the test function in the
volume terms:

| Operator on the test function | Boundary traces                                        |
|-------------------------------|--------------------------------------------------------|
| value ``δu``                  | none                                                   |
| gradient ``∇δu``              | ``δu``                                                 |
| second gradient ``∇∇δu``      | ``δu``, ``∂_{n_0} δu`` on faces, ``δu`` on sharp edges |
| divergence ``Div δq``         | ``δq \cdot n_0``                                       |

So the admissible traces of each field are computed, not declared. This also gets the
roles right in mixed formulations: for mixed Darcy flow the normal flux is the essential
condition and the pressure the natural one, simply because the flux test function appears
under a divergence.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------- Traces of test operators ------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    boundary_traces(vertex)

Return the boundary traces left by integrating by parts a volume term whose test side is
the variation of `vertex`.
"""
boundary_traces(::AbstractVertex) = ()
boundary_traces(::Value) = ()
boundary_traces(::ReferenceGradient) = (ValueTrace(),)
boundary_traces(::Divergence) = (NormalComponentTrace(),)
boundary_traces(::SecondGradient) = (ValueTrace(), NormalGradientTrace(), EdgeValueTrace())
boundary_traces(::AbstractKinematicVertex) = (ValueTrace(),)
boundary_traces(vertex::Rate) = boundary_traces(vertex.vertex)
boundary_traces(vertex::PullBack) = boundary_traces(vertex.vertex)
boundary_traces(vertex::PushForward) = boundary_traces(vertex.vertex)


"""
    derived_traces(problem, field)

Return the boundary traces of the field called `field`, derived from the operators that
act on its test function in the volume terms of `problem`.
"""
function derived_traces(problem::Problem, field::Symbol)
    traces = AbstractTrace[]
    for (equation, term) in problem_terms(problem)
        # Only volume terms tested by this field are integrated by parts.
        if test_field(equation) === field && is_volume_term(term)
            for trace in boundary_traces(term_test(term).vertex)
                push_unique!(traces, trace)
            end
        end
    end
    return traces
end


# ---------------------------------------------------------------------------------------- #
# -------------------------------- Status of every trace --------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractTraceStatus

Supertype of the statuses a derived trace can have on a boundary.
"""
abstract type AbstractTraceStatus end


"""
    Prescribed()

The trace is prescribed by an essential boundary condition.
"""
struct Prescribed <: AbstractTraceStatus end


"""
    Loaded()

The conjugate boundary quantity is prescribed by a natural boundary condition.
"""
struct Loaded <: AbstractTraceStatus end


"""
    DefaultedToZero()

Nothing is prescribed, so the homogeneous natural condition applies (zero traction, zero
double force, …).
"""
struct DefaultedToZero <: AbstractTraceStatus end


"""
    BoundaryTraceStatus(boundary, field, trace, status)

The status of one derived trace of one field on one boundary.

# Fields
- `boundary::Symbol`: the name of the facet set.
- `field::Symbol`: the field.
- `trace::T`: the derived [`AbstractTrace`](@ref).
- `status::S`: the [`AbstractTraceStatus`](@ref): whether the trace is prescribed, loaded
  or defaulted.
"""
struct BoundaryTraceStatus{T<:AbstractTrace,S<:AbstractTraceStatus}
    boundary::Symbol
    field::Symbol
    trace::T
    status::S
end


"""
    boundary_status(problem)

Return the status of every derived trace of every field on every boundary that `problem`
names in a boundary condition.

Only boundaries named in the problem are known before the mesh is; a complete check over
every facet set of the mesh belongs to the compile step.
"""
function boundary_status(problem::Problem)
    statuses = BoundaryTraceStatus[]

    # We collect the boundaries in the order the conditions first name them.
    boundaries = Symbol[]
    for bc in problem.boundary_conditions
        push_unique!(boundaries, bc_boundary(bc))
    end

    # For every boundary, field and derived trace we look for a condition that covers it.
    for boundary in boundaries
        for field in field_names(problem)
            for trace in derived_traces(problem, field)
                status = trace_status(problem, boundary, field, trace)
                push!(statuses, BoundaryTraceStatus(boundary, field, trace, status))
            end
        end
    end
    return statuses
end


"""
    trace_status(problem, boundary, field, trace)

Return whether `trace` of `field` on `boundary` is prescribed, loaded or defaulted to zero,
after translating every boundary condition through the problem's template.
"""
function trace_status(
    problem::Problem,
    boundary::Symbol,
    field::Symbol,
    trace::AbstractTrace,
)
    for bc in problem.boundary_conditions
        # Conditions written on the continuum theory are mapped onto the discretised fields.
        target_field, target_trace =
            translate_trace(problem.template, bc_field(bc), bc_trace(bc))
        if bc_boundary(bc) === boundary && target_field === field && target_trace == trace
            if bc isa EssentialBC
                return Prescribed()
            else
                return Loaded()
            end
        end
    end
    return DefaultedToZero()
end
