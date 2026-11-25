#=
Quadrature.

The quadrature order of a cell is chosen from the interpolations (and the nonlinearity of
the integrand) of the fields it carries, as a function of the interpolation rather than of
a bare integer, so that families other than Lagrange choose correctly when they are
admitted.
=#


"""
    quadrature_order(interpolations, problem)

Return the quadrature order for cells carrying `interpolations` in `problem`.
"""
function quadrature_order(interpolations, problem::Problem)
    not_yet_implemented("selection of the quadrature order")
    return nothing
end
