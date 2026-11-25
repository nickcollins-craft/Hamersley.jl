#=
Output.

Results are written by evaluating (or projecting) fields and quadrature-point quantities
at output points, never by reshaping the solution vector node by node. That keeps output
correct for interpolations whose degrees of freedom are not nodal values.
=#


"""
    write_output(filename, compiled, solution, state)

Write the fields of `solution` and the quadrature-point quantities in `state` for
`compiled` to `filename` in VTK format.
"""
function write_output(filename::AbstractString, compiled::CompiledProblem, solution, state)
    not_yet_implemented("writing output")
    return nothing
end
