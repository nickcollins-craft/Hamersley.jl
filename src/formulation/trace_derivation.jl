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
    for (equation, term) in weak_form_terms(problem)
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
    PrescribedTrace()

The trace is prescribed by an essential boundary condition.
"""
struct PrescribedTrace <: AbstractTraceStatus end


"""
    LoadedTrace()

The conjugate boundary quantity is prescribed by a natural boundary condition.
"""
struct LoadedTrace <: AbstractTraceStatus end


"""
    TermLoadedTrace()

The trace is loaded by a boundary term of the weak form: a Robin condition, or a boundary
term written in an equation (a convective face, say).
"""
struct TermLoadedTrace <: AbstractTraceStatus end


"""
    DefaultedToZero()

Nothing is prescribed, so the homogeneous natural condition applies (zero traction, zero
double force, …).
"""
struct DefaultedToZero <: AbstractTraceStatus end


"""
    status_label(status)

Return the word that reports print for `status`.
"""
status_label(::PrescribedTrace) = "Prescribed"
status_label(::LoadedTrace) = "Loaded"
status_label(::TermLoadedTrace) = "LoadedByTerm"
status_label(::DefaultedToZero) = "DefaultedToZero"


"""
    BoundaryTraceStatus(region, field, trace, components, status)

The status of some components of one derived trace of one field on one region.

# Fields
- `region::R`: the [`AbstractRegion`](@ref).
- `field::Symbol`: the field.
- `trace::T`: the derived [`AbstractTrace`](@ref).
- `components::C`: the [`AbstractComponentSelection`](@ref) the status applies to.
- `status::S`: the [`AbstractTraceStatus`](@ref): whether those components are
  prescribed, loaded or defaulted.
"""
struct BoundaryTraceStatus{
    R<:AbstractRegion,
    T<:AbstractTrace,
    C<:AbstractComponentSelection,
    S<:AbstractTraceStatus,
}
    region::R
    field::Symbol
    trace::T
    components::C
    status::S
end


"""
    condition_status(bc)

Return the status that `bc` gives the trace it acts on.
"""
condition_status(::EssentialBC) = PrescribedTrace()
condition_status(::NaturalBC) = LoadedTrace()
condition_status(::RobinBC) = TermLoadedTrace()


"""
    TraceClaim(region, field, trace, components, status, source)

One claim on some components of a trace on a region, by a boundary condition or by a
boundary term written in an equation, after the template's translation. Trace statuses
and the conflict rule are both built from these.

# Fields
- `region::AbstractRegion`: the region.
- `field::Symbol`: the (translated) field.
- `trace::AbstractTrace`: the (translated) trace.
- `components::AbstractComponentSelection`: the components claimed.
- `status::AbstractTraceStatus`: the status the claim gives.
- `source::Any`: the condition or term making the claim.
"""
struct TraceClaim
    region::AbstractRegion
    field::Symbol
    trace::AbstractTrace
    components::AbstractComponentSelection
    status::AbstractTraceStatus
    source::Any
end


"""
    trace_claims(problem, conditions)

Return the claims that `conditions` and the boundary terms written in the equations of
`problem` make on the traces of the discretised fields. Discrete conditions (a pinned
vertex, say) make no claim on any trace of the weak form, but with
`include_discrete = true` they are included, so that two conditions on the same points can
be found to conflict.
"""
function trace_claims(problem::Problem, conditions; include_discrete::Bool=false)
    claims = TraceClaim[]

    # Every condition claims its translated trace, on its components.
    for bc in conditions
        if is_discrete_condition(bc) && !include_discrete
            continue
        end
        field, trace = translate_trace(problem.template, bc_field(bc), bc_trace(bc))
        claim = TraceClaim(
            bc_region(bc),
            field,
            trace,
            bc_components(bc),
            condition_status(bc),
            bc,
        )
        push!(claims, claim)
    end

    # A boundary term restricted to a facet set, varying the value of its equation's field,
    # loads that field's value trace there.
    for (equation, term) in weak_form_terms(problem)
        domain = term_domain(term)
        if is_boundary_term(term) && domain !== nothing && term_test(term).vertex isa Value
            claim = TraceClaim(
                FacetSet(domain),
                test_field(equation),
                ValueTrace(),
                AllComponents(),
                TermLoadedTrace(),
                term,
            )
            push!(claims, claim)
        end
    end
    return claims
end


"""
    whole_boundary_term_fields(problem)

Return the fields whose value trace is loaded on every boundary by a boundary term that is
not restricted to a facet set.
"""
function whole_boundary_term_fields(problem::Problem)
    fields = Symbol[]
    for (equation, term) in weak_form_terms(problem)
        if is_boundary_term(term) &&
           term_domain(term) === nothing &&
           term_test(term).vertex isa Value
            push_unique!(fields, test_field(equation))
        end
    end
    return fields
end


"""
    boundary_status(problem, phase = 1)

Return the status of every derived trace of every field on every facet or edge region
that `problem` names (in a condition of the phase, or as the domain of a boundary term),
component by component, in the phase `phase` (a name or a position).

Which components remain once some are claimed is only known with the spatial dimension,
so the remainder is reported symbolically (for example "all components except 3"); a
complete check over every facet set of the mesh belongs to the compile step.
"""
function boundary_status(problem::Problem, phase=1)
    statuses = BoundaryTraceStatus[]
    claims = trace_claims(problem, phase_conditions(problem, phase))
    everywhere = whole_boundary_term_fields(problem)

    # We collect the regions in the order the claims first name them.
    regions = AbstractRegion[]
    for claim in claims
        push_unique!(regions, claim.region)
    end

    # For every region, field and derived trace we report each claim, then what is left.
    for region in regions
        for field in field_names(problem)
            for trace in derived_traces(problem, field)
                if !carries_trace(region, trace)
                    continue
                end
                append!(statuses, trace_statuses(claims, everywhere, region, field, trace))
            end
        end
    end
    return statuses
end


"""
    trace_statuses(claims, everywhere, region, field, trace)

Return the statuses of `trace` of `field` on `region`: one per claim on it, and one for the
components that no claim covers, which default to zero (unless a boundary term on the
whole boundary loads the value trace of `field`).
"""
function trace_statuses(claims, everywhere, region, field, trace)
    statuses = BoundaryTraceStatus[]
    remaining = AllComponents()
    for claim in claims
        if claim.region == region && claim.field === field && claim.trace == trace
            push!(
                statuses,
                BoundaryTraceStatus(region, field, trace, claim.components, claim.status),
            )
            remaining = remaining === nothing ? nothing :
                remove_components(remaining, claim.components)
        end
    end

    # What no claim covers is loaded by a whole-boundary term, or defaults to zero.
    if remaining !== nothing
        status = if trace isa ValueTrace && field in everywhere
            TermLoadedTrace()
        else
            DefaultedToZero()
        end
        push!(statuses, BoundaryTraceStatus(region, field, trace, remaining, status))
    end
    return statuses
end
