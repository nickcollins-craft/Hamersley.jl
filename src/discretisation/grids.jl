#=
Grids.

Hamersley takes its meshes as Ferrite grids with named cell, facet and vertex sets,
wherever they come from: the `Meshing` submodule, Ferrite's own `generate_grid`, or a
user's file. This file checks that a grid provides every set and region a problem refers
to.
=#


"""
    grid_issues(grid, problem)

Return a message for every cell set, facet set or vertex set that `problem` names but
`grid` does not provide, and for every edge region that `grid` cannot provide.
"""
function grid_issues(grid::Ferrite.AbstractGrid, problem::Problem)
    messages = String[]

    # Boundary conditions name regions, in every phase (releases included).
    for phase in problem_phases(problem)
        for entry in phase.conditions
            append!(messages, region_grid_issues(grid, bc_region(entry)))
        end
    end

    # Restricted terms name cell sets (volume terms) or facet sets (boundary terms).
    facet_set_names = keys(Ferrite.getfacetsets(grid))
    cell_set_names = keys(Ferrite.getcellsets(grid))
    for (_, term) in weak_form_terms(problem)
        domain = term_domain(term)
        if domain === nothing
            continue
        elseif is_volume_term(term) && !(String(domain) in cell_set_names)
            push!(messages, "the grid has no cell set named $domain")
        elseif is_boundary_term(term) && !(String(domain) in facet_set_names)
            push!(messages, "the grid has no facet set named $domain")
        end
    end

    # A region named in several phases is reported once.
    return unique(messages)
end


"""
    region_grid_issues(grid, region)

Return a message if `grid` does not provide `region`. Ferrite keys its sets by strings,
while Hamersley names them with symbols.
"""
function region_grid_issues(grid::Ferrite.AbstractGrid, region::FacetSet)
    messages = String[]
    if !(String(region.name) in keys(Ferrite.getfacetsets(grid)))
        push!(messages, "the grid has no facet set named $(region.name)")
    end
    return messages
end

function region_grid_issues(grid::Ferrite.AbstractGrid, region::VertexSet)
    messages = String[]
    if !(String(region.name) in keys(Ferrite.getvertexsets(grid)))
        push!(messages, "the grid has no vertex set named $(region.name)")
    end
    return messages
end

function region_grid_issues(grid::Ferrite.AbstractGrid, region::EdgeSet)
    messages = String[]

    # Edges are a region of their own only in 3D; in 2D the facets are the edges.
    if Ferrite.getspatialdim(grid) != 3
        push!(
            messages,
            "$(region_label(region)) needs a 3D grid; in 2D the facets are the edges, so " *
                "use a FacetSet",
        )
        return messages
    end

    # Ferrite grids carry no named edge sets, so only edges between facet sets can be
    # found on one.
    if region.between === nothing
        push!(
            messages,
            "the grid cannot provide $(region_label(region)): Ferrite grids carry no " *
                "named edge sets, so give the edge as the meeting of two facet sets, " *
                "EdgeSet(:first, :second)",
        )
    else
        for name in region.between
            append!(messages, region_grid_issues(grid, FacetSet(name)))
        end
    end
    return messages
end


"""
    validate_grid(grid, problem)

Throw an `ArgumentError` listing every set that `problem` names but `grid` lacks.
"""
function validate_grid(grid::Ferrite.AbstractGrid, problem::Problem)
    # We report every missing set at once, rather than only the first.
    messages = grid_issues(grid, problem)
    if !isempty(messages)
        throw(ArgumentError(join(messages, "\n")))
    end
    return nothing
end
