#=
Regions and component selections.

A boundary condition acts on a region of the boundary (a set of facets, edges or
vertices) and, for a trace with several components, on a selection of those components.
Regions are named by the sets of the mesh, but the formulation only needs their kind:
facet regions carry the traces that integration by parts produces, edge regions carry the
edge traces of second-gradient continua, and vertex regions only ever carry discrete
conditions (a pinned point, say).

Which components exist depends on the spatial dimension, which is only known once the mesh
is. A selection is therefore stored symbolically, as every component, some of them, or
every component except some of them, and the set operations below never need the
dimension.
=#


# ---------------------------------------------------------------------------------------- #
# --------------------------------------- Regions ---------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractRegion

Supertype of the regions of the boundary that a condition can act on.
"""
abstract type AbstractRegion end


"""
    FacetSet(name)

The facets (faces in 3D, edges in 2D) of the mesh's facet set `name`. A bare `Symbol`
given as a region means a facet set.

# Fields
- `name::Symbol`: the name of the facet set.
"""
struct FacetSet <: AbstractRegion
    name::Symbol
end


"""
    EdgeSet(name)
    EdgeSet(first_facet_set, second_facet_set)

The edges of a 3D mesh: either a named edge set, or the edges where two facet sets meet
(for example `EdgeSet(:top, :front)`). Ferrite grids carry no named edge sets, so the
second form is the one that can be checked against a grid today.

# Fields
- `name::Union{Symbol,Nothing}`: the name of the edge set, or `nothing`.
- `between::Union{Tuple{Symbol,Symbol},Nothing}`: the two facet sets whose shared edges
  form the region, or `nothing`.
"""
struct EdgeSet <: AbstractRegion
    name::Union{Symbol,Nothing}
    between::Union{Tuple{Symbol,Symbol},Nothing}
end

EdgeSet(name::Symbol) = EdgeSet(name, nothing)
EdgeSet(first_facet_set::Symbol, second_facet_set::Symbol) =
    EdgeSet(nothing, (first_facet_set, second_facet_set))


"""
    VertexSet(name)

The vertices of the mesh's vertex set `name`, for conditions at individual points (most
often pinning a point to remove a rigid-body mode).

# Fields
- `name::Symbol`: the name of the vertex set.
"""
struct VertexSet <: AbstractRegion
    name::Symbol
end


"""
    as_region(on)

Return the region that `on` names: a region is returned as it is, and a `Symbol` is a
[`FacetSet`](@ref).
"""
as_region(region::AbstractRegion) = region
as_region(name::Symbol) = FacetSet(name)


"""
    region_label(region)

Return a short description of `region` for messages and reports. A facet set is labelled
by its name alone, as it always has been.
"""
region_label(region::FacetSet) = string(region.name)
region_label(region::VertexSet) = "vertex set $(region.name)"

function region_label(region::EdgeSet)
    # A named edge set is labelled by its name, otherwise by the two facet sets it joins.
    if region.name !== nothing
        return "edge set $(region.name)"
    else
        first_facet_set, second_facet_set = region.between
        return "edge between $first_facet_set and $second_facet_set"
    end
end


"""
    facet_set_name(region)

Return the name of the facet set that `region` is, or `nothing` for any other region. A
weak term restricted to a domain names its facet set this way.
"""
facet_set_name(region::FacetSet) = region.name
facet_set_name(::AbstractRegion) = nothing


# ---------------------------------------------------------------------------------------- #
# --------------------------------- Component selections --------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AbstractComponentSelection

Supertype of the selections of the components of a trace (or of a boundary quantity) that
a condition acts on.
"""
abstract type AbstractComponentSelection end


"""
    AllComponents()

Every component: the default.
"""
struct AllComponents <: AbstractComponentSelection end


"""
    SomeComponents(indices)

The components `indices`, numbered from one, sorted and without repetition.

# Fields
- `indices::Tuple{Vararg{Int}}`: the selected components.
"""
struct SomeComponents <: AbstractComponentSelection
    indices::Tuple{Vararg{Int}}
end


"""
    AllComponentsExcept(indices)

Every component except `indices`. This is what remains of a condition on every component
once a later phase replaces some of them, which can be stated without knowing how many
components there are.

# Fields
- `indices::Tuple{Vararg{Int}}`: the excluded components.
"""
struct AllComponentsExcept <: AbstractComponentSelection
    indices::Tuple{Vararg{Int}}
end


"""
    component_selection(components)

Return the selection that the `components` keyword of a condition names: `nothing` (or
`AllComponents()`) for every component, an integer for one, or a collection of integers.
Indices are sorted and repeats are dropped; their validity is checked by the formulation
verification (see [`selection_issues`](@ref)).
"""
component_selection(::Nothing) = AllComponents()
component_selection(selection::AbstractComponentSelection) = selection
component_selection(index::Integer) = SomeComponents((Int(index),))
component_selection(indices) = SomeComponents(Tuple(sort!(unique!(collect(Int, indices)))))


"""
    selection_issues(selection)

Return a message for every way in which `selection` is malformed (no components, or
indices below one).
"""
selection_issues(::AbstractComponentSelection) = String[]

function selection_issues(selection::Union{SomeComponents,AllComponentsExcept})
    messages = String[]

    # A selection of nothing selects nothing, which is never what was meant.
    if selection isa SomeComponents && isempty(selection.indices)
        push!(messages, "selects no components")
    end

    # Components are numbered from one.
    if any(<(1), selection.indices)
        push!(messages, "selects components numbered below one")
    end
    return messages
end


"""
    selects_all(selection)

Return `true` when `selection` is every component.
"""
selects_all(::AbstractComponentSelection) = false
selects_all(::AllComponents) = true


"""
    selected_count(selection)

Return the number of selected components when it is known without the spatial dimension,
or `nothing`.
"""
selected_count(::AbstractComponentSelection) = nothing
selected_count(selection::SomeComponents) = length(selection.indices)


"""
    overlaps(a, b)

Return `true` when the selections `a` and `b` may share a component.

Two selections that each exclude finitely many components always share one in a space of
high enough dimension, so they are taken to overlap; a selection of some components
overlaps one that excludes some exactly when it selects a component that is not excluded.
"""
overlaps(::AbstractComponentSelection, ::AbstractComponentSelection) = true
overlaps(a::SomeComponents, b::SomeComponents) = !isempty(intersect(a.indices, b.indices))
overlaps(a::SomeComponents, b::AllComponentsExcept) = !issubset(a.indices, b.indices)
overlaps(a::AllComponentsExcept, b::SomeComponents) = overlaps(b, a)


"""
    remove_components(selection, removed)

Return the selection that is left of `selection` once the components of `removed` are
taken away, or `nothing` when nothing is left. This is how a condition inherited from an
earlier phase is restricted when a later phase replaces some of its components.
"""
function remove_components(selection::AbstractComponentSelection, removed::AllComponents)
    # Removing every component leaves nothing, whatever was selected.
    return nothing
end

function remove_components(::AllComponents, removed::SomeComponents)
    remaining = AllComponentsExcept(removed.indices)
    return remaining
end

function remove_components(selection::AllComponentsExcept, removed::SomeComponents)
    # The excluded components grow by the removed ones.
    excluded = Tuple(sort!(union(collect(selection.indices), collect(removed.indices))))
    remaining = AllComponentsExcept(excluded)
    return remaining
end

function remove_components(selection::SomeComponents, removed::SomeComponents)
    # Only the selected components that are not removed remain.
    kept = Tuple(setdiff(selection.indices, removed.indices))
    if isempty(kept)
        return nothing
    else
        return SomeComponents(kept)
    end
end

function remove_components(selection::SomeComponents, removed::AllComponentsExcept)
    # Only the selected components that the removal leaves alone remain.
    kept = Tuple(intersect(selection.indices, removed.indices))
    if isempty(kept)
        return nothing
    else
        return SomeComponents(kept)
    end
end

function remove_components(selection::AllComponents, removed::AllComponentsExcept)
    # Removing every component but A from every component leaves exactly A.
    return remaining_selection(removed.indices)
end

function remove_components(selection::AllComponentsExcept, removed::AllComponentsExcept)
    # Removing every component but A from every component but B leaves A \ B.
    kept = Tuple(setdiff(removed.indices, selection.indices))
    return remaining_selection(kept)
end


"""
    remaining_selection(indices)

Return `SomeComponents(indices)`, or `nothing` when `indices` is empty.
"""
function remaining_selection(indices::Tuple)
    if isempty(indices)
        return nothing
    else
        return SomeComponents(indices)
    end
end


"""
    selection_label(selection)

Return a short description of `selection` for messages and reports, or an empty string
for every component.
"""
selection_label(::AllComponents) = ""
selection_label(selection::SomeComponents) =
    " [components $(join(string.(selection.indices), ", "))]"
selection_label(selection::AllComponentsExcept) =
    " [all components except $(join(string.(selection.indices), ", "))]"
