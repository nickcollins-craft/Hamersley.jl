#=
Assembly modes and the element kernel contract.

Every element kernel can be asked for the residual only, the tangent only, or both. The
residual-only mode is what makes linear problems cheap after the first step: their
tangent is assembled and factorised once, and every later step assembles only the
right-hand side.

An element kernel follows three rules so that elements other than Lagrange can be
admitted later without rewriting it:

1. It evaluates fields through Ferrite's cell values (`shape_value`, `shape_gradient`,
   `function_value`, …), never by treating degrees of freedom as nodal values, and never
   assumes that the number of basis functions equals nodes times components.
2. It reinitialises cell values with the cell itself (a `CellCache`), not only with its
   coordinates, because ``H(\mathrm{div})`` and ``H(\mathrm{curl})`` elements need the
   cell's orientation.
3. It loops over the integration domains its terms name (cells, boundary facets, and
   later interior facets) rather than hard-coding cells and boundaries.
=#


"""
    AbstractAssemblyMode

Supertype of the assembly modes.
"""
abstract type AbstractAssemblyMode end


"""
    ResidualOnly()

Assemble the residual, but not the tangent.
"""
struct ResidualOnly <: AbstractAssemblyMode end


"""
    TangentOnly()

Assemble the tangent, but not the residual.
"""
struct TangentOnly <: AbstractAssemblyMode end


"""
    ResidualAndTangent()

Assemble both the residual and the tangent.
"""
struct ResidualAndTangent <: AbstractAssemblyMode end


needs_residual(::ResidualOnly) = true
needs_residual(::TangentOnly) = false
needs_residual(::ResidualAndTangent) = true

needs_tangent(::ResidualOnly) = false
needs_tangent(::TangentOnly) = true
needs_tangent(::ResidualAndTangent) = true


"""
    element_kernel!(local_residual, local_tangent, mode, kernel, cell, local_solution,
                    state)

Evaluate the element contributions of `kernel` on `cell` in the given `mode`, overwriting
`local_residual` and/or `local_tangent` (whichever `mode` asks for) and updating the
internal state of the cell's quadrature points in `state`.
"""
function element_kernel!(
    local_residual,
    local_tangent,
    mode::AbstractAssemblyMode,
    kernel,
    cell,
    local_solution,
    state,
)
    not_yet_implemented("element kernels")
    return nothing
end
