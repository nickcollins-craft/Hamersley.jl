#=
Grids.

Hamersley takes its meshes as Ferrite grids with named cell and facet sets, wherever they
come from: the `Meshing` submodule, Ferrite's own `generate_grid`, or a user's file. This
file checks that a grid provides every set a problem refers to.
=#


"""
    grid_issues(grid, problem)

Return a message for every cell set or facet set that `problem` names but `grid` does not
provide.
"""
function grid_issues(grid::Ferrite.AbstractGrid, problem::Problem)
    messages = String[]

    # Ferrite keys its sets by strings, while Hamersley names them with symbols.
    facet_set_names = keys(Ferrite.getfacetsets(grid))
    cell_set_names = keys(Ferrite.getcellsets(grid))

    # Boundary conditions name facet sets.
    for bc in problem.boundary_conditions
        if !(String(bc_boundary(bc)) in facet_set_names)
            push!(messages, "the grid has no facet set named $(bc_boundary(bc))")
        end
    end

    # Restricted terms name cell sets (volume terms) or facet sets (boundary terms).
    for (_, term) in problem_terms(problem)
        domain = term_domain(term)
        if domain === nothing
            continue
        elseif is_volume_term(term) && !(String(domain) in cell_set_names)
            push!(messages, "the grid has no cell set named $domain")
        elseif is_boundary_term(term) && !(String(domain) in facet_set_names)
            push!(messages, "the grid has no facet set named $domain")
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
