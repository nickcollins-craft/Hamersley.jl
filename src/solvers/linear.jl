#=
The cached-factorisation strategy for linear problems.

For a problem whose residual is affine, one Newton step from the predictor is exact:

```math
x = x_{\text{pred}} - K^{-1} R(x_{\text{pred}}).
```

So the tangent ``K`` is assembled, constrained and factorised once, and reused for as long
as its `FactorisationKey` (time step, parameters, constrained set) is unchanged. Every
step then assembles the residual only, applies the Dirichlet values to the right-hand side
with Ferrite's `get_rhs_data` / `apply_rhs!` (which leave the factorised matrix alone),
and back-substitutes. No material state is updated, no convergence check is made.
=#


"""
    CachedFactorisationStrategy()

Solve an affine problem with a single back-substitution per step, refactorising the
tangent only when its [`FactorisationKey`](@ref) changes.
"""
struct CachedFactorisationStrategy <: AbstractSolverStrategy end
