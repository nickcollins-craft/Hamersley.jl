#=
Global assembly.

Global assembly loops over the integration domains of a compiled problem, calls the
element kernels in the requested mode, and scatters their contributions into the global
residual and tangent. Only the blocks the solver strategy asks for are assembled, so that
staggered schemes can skip off-diagonal blocks of a multiphysics tangent.
=#


"""
    assemble!(residual, tangent, mode, compiled, solution, state)

Assemble the global `residual` and/or `tangent` (whichever `mode` asks for) of `compiled`
at `solution`, updating the internal state in `state`. The arrays are overwritten.
"""
function assemble!(
    residual,
    tangent,
    mode::AbstractAssemblyMode,
    compiled::CompiledProblem,
    solution,
    state,
)
    not_yet_implemented("global assembly")
    return nothing
end
