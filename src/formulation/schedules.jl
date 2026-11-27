#=
Schedules.

Boundary conditions may change during an analysis: a face loaded by a traction until some
time and then held where it stopped, say. A schedule declares this as a sequence of
phases, each ending with an explicit trigger, rather than hiding it in callbacks, so that
every phase can be verified as a problem in its own right before the analysis starts.

Each phase inherits the conditions of the previous one. A condition in a phase replaces
the inherited conditions on the same key (region, field, trace), but only on the
components it shares with them: an essential condition on ``u_z`` replaces the ``z``
component of an inherited traction and leaves the others loaded. `Release` returns a key
to the homogeneous natural default.

A problem with a plain tuple of boundary conditions has a single phase that never ends.
=#


# ---------------------------------------------------------------------------------------- #
# --------------------------------------- Triggers --------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractTrigger

Supertype of the events that end a phase.
"""
abstract type AbstractTrigger end


"""
    AtTime(time)

End the phase at the time `time`. The step controller shortens the last step of the phase
so that it lands exactly on it.

# Fields
- `time::Float64`: the time at which the phase ends.
"""
struct AtTime <: AbstractTrigger
    time::Float64
end


"""
    Never()

The phase never ends: it is the last one, and runs to the end of the analysis.
"""
struct Never <: AbstractTrigger end


# ---------------------------------------------------------------------------------------- #
# ---------------------------------------- Phases ---------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    Release(field; on, trace = ValueTrace(), components = nothing)

Return the trace `trace` of `field` on the region `on` (or some of its components) to the
homogeneous natural default, removing whatever condition a phase inherited there.

# Fields
- `field::Symbol`: the field.
- `trace::T`: the [`AbstractTrace`](@ref).
- `region::R`: the [`AbstractRegion`](@ref).
- `components::C`: the released components.
"""
struct Release{T<:AbstractTrace,R<:AbstractRegion,C<:AbstractComponentSelection}
    field::Symbol
    trace::T
    region::R
    components::C
end


function Release(
    field::Symbol;
    on::Union{Symbol,AbstractRegion},
    trace::AbstractTrace=ValueTrace(),
    components=nothing,
)
    release = Release(field, trace, as_region(on), component_selection(components))
    return release
end


bc_field(release::Release) = release.field
bc_trace(release::Release) = release.trace
bc_region(release::Release) = release.region
bc_components(release::Release) = release.components


"""
    Phase(name; conditions = (), until = Never())

A phase of an analysis: the boundary conditions that change at its start (new conditions
and [`Release`](@ref)s) and the trigger that ends it.

# Fields
- `name::Symbol`: a name for the phase, used in messages and reports.
- `conditions::Cs`: a tuple of boundary conditions and releases.
- `until::T`: the [`AbstractTrigger`](@ref) that ends the phase.
"""
struct Phase{Cs<:Tuple,T<:AbstractTrigger}
    name::Symbol
    conditions::Cs
    until::T
end


Phase(name::Symbol; conditions::Tuple=(), until::AbstractTrigger=Never()) =
    Phase(name, conditions, until)

phase_name(phase::Phase) = phase.name


"""
    Schedule(phases...)

A sequence of [`Phase`](@ref)s.

# Fields
- `phases::Ps`: the phases, in order.
"""
struct Schedule{Ps<:Tuple}
    phases::Ps
end


Schedule(phases::Phase...) = Schedule(phases)


"""
    single_phase_schedule(conditions)

Return the schedule of a problem given a plain tuple of boundary conditions: one phase,
called `default`, that never ends.
"""
single_phase_schedule(conditions::Tuple) = Schedule(Phase(:default; conditions=conditions))


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- Inheritance ------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    same_key(a, b)

Return `true` when the conditions (or releases) `a` and `b` act on the same field, trace
and region, so that they compete for the same components.
"""
function same_key(a, b)
    same_field = bc_field(a) === bc_field(b)
    same_trace = bc_trace(a) == bc_trace(b)
    same_region = bc_region(a) == bc_region(b)
    return same_field && same_trace && same_region
end


"""
    restrict_inherited(inherited, entry)

Return the conditions of `inherited` with the components of `entry` (a new condition or a
release on the same key) taken away; a condition left with no components is dropped. The
second return value says whether `entry` replaced or released anything.
"""
function restrict_inherited(inherited::Vector{AbstractBoundaryCondition}, entry)
    kept = AbstractBoundaryCondition[]
    touched = false
    for bc in inherited
        # Conditions on other keys, or on other components, are left as they are.
        if !same_key(bc, entry) || !overlaps(bc_components(bc), bc_components(entry))
            push!(kept, bc)
            continue
        end

        # Otherwise only the components that the entry does not take over remain.
        touched = true
        remaining = remove_components(bc_components(bc), bc_components(entry))
        if remaining !== nothing
            push!(kept, with_components(bc, remaining))
        end
    end
    return kept, touched
end


"""
    resolved_phases(schedule)

Return, for every phase of `schedule`, the vector of boundary conditions that apply in it
once inheritance, replacement and releases are taken into account.
"""
function resolved_phases(schedule::Schedule)
    resolved = Vector{AbstractBoundaryCondition}[]
    inherited = AbstractBoundaryCondition[]
    for phase in schedule.phases
        # The phase's entries replace or release inherited conditions on their keys.
        kept = inherited
        for entry in phase.conditions
            kept, _ = restrict_inherited(kept, entry)
        end

        # The phase's own conditions are added to what is kept, and passed on.
        own = AbstractBoundaryCondition[
            entry for entry in phase.conditions if entry isa AbstractBoundaryCondition
        ]
        conditions = vcat(kept, own)
        push!(resolved, conditions)
        inherited = conditions
    end
    return resolved
end


"""
    ineffective_releases(schedule)

Return, as `(phase, release)` pairs, every release that removes nothing it inherited
(usually a mistyped region or field).
"""
function ineffective_releases(schedule::Schedule)
    pairs = Tuple{Phase,Release}[]
    inherited = AbstractBoundaryCondition[]
    resolved = resolved_phases(schedule)
    for (index, phase) in enumerate(schedule.phases)
        # A release is checked against what the phase inherited, not against its own
        # conditions.
        for entry in phase.conditions
            if entry isa Release
                _, touched = restrict_inherited(inherited, entry)
                if !touched
                    push!(pairs, (phase, entry))
                end
            end
        end
        inherited = resolved[index]
    end
    return pairs
end
