#=
Degrees of freedom.

The degree-of-freedom handler is built per subdomain (one Ferrite `SubDofHandler` per cell
set) from the start, even when there is only one subdomain. That is what will later allow
mixed element shapes in one mesh, and fields that only exist in part of the domain (a pore
pressure in one layer, say).
=#


"""
    build_dofhandler(problem, grid)

Build the Ferrite `DofHandler` for `problem` on `grid`, with one `SubDofHandler` per
subdomain and the fields added in the order `problem` declares them.
"""
function build_dofhandler(problem::Problem, grid::Ferrite.AbstractGrid)
    not_yet_implemented("construction of the degree-of-freedom handler")
    return nothing
end
