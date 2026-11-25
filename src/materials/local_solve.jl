#=
Local (material-point) nonlinear solves.

Return-mapping algorithms solve a small nonlinear system at every quadrature point, many
millions of times per global iteration. We hide the choice of local solver behind the
single function `solve_local`, so that it can be changed without touching the material
library.

We deliberately do not depend on NonlinearSolve here. When the first plasticity models are
written we will compare a hand-written Newton iteration on static arrays (allocation-free,
returning the local Jacobian needed for the consistent tangent) against
SimpleNonlinearSolve, and choose then.
=#


"""
    AbstractLocalSolveStatus

Supertype of the outcomes of a local solve. A failed local solve must not throw: it
signals the global solver to cut the step back.
"""
abstract type AbstractLocalSolveStatus end


"""
    LocalSolveConverged()

The local solve met its tolerances.
"""
struct LocalSolveConverged <: AbstractLocalSolveStatus end


"""
    LocalSolveFailed()

The local solve did not converge; the global step should be cut back.
"""
struct LocalSolveFailed <: AbstractLocalSolveStatus end


"""
    solve_local(residual, x_nought, parameters; relative_tolerance, absolute_tolerance)

Solve the local system ``R(x; \\text{parameters}) = 0`` starting from `x_nought`.

# Returns
A tuple `(x, status, ∂R_∂x)` with the converged solution, an
[`AbstractLocalSolveStatus`](@ref), and the Jacobian ``∂R/∂x`` at the solution. The
Jacobian is what the implicit function theorem needs to form the consistent tangent
``\\mathrm{d}x/\\mathrm{d}(\\text{inputs}) = -(∂R/∂x)^{-1} ∂R/∂(\\text{inputs})``.
"""
function solve_local(
    residual,
    x_nought,
    parameters;
    relative_tolerance::Real = 1e-12,
    absolute_tolerance::Real = 1e-12,
)
    not_yet_implemented("the material-point nonlinear solver (solve_local)")
    return nothing
end
