#=
The cached-factorisation strategy for linear problems.

For a problem whose residual is affine, one Newton step from the predictor is exact:

```math
x = x_{\text{pred}} - K^{-1} R(x_{\text{pred}}).
```

So the tangent ``K`` is assembled, constrained and factorised once, and reused for as long
as its `FactorisationKey` (time step, parameters, constrained set) is unchanged. Every
step then:

1. assembles the residual only at the predictor, with the laws reading the converged state
   of the previous step;
2. applies the essential values to the right-hand side with Ferrite's `get_rhs_data` /
   `apply_rhs!` (which leave the factorised matrix alone), and back-substitutes;
3. when the problem has internal state (`has_internal_state`), assembles the residual once
   more at the solution, which is what writes the new state at every quadrature point,
   and commits it.

Being affine does not make a problem stateless: a linear Maxwell model needs its stress
``σ_n`` carried from step to step, and skipping step 3 would freeze it at its initial value,
so that the material would never relax. The residual of step 3 should vanish to round-off,
so its norm is a free check that the affine-ness the laws declare is true. No convergence
loop is run.
=#


"""
    CachedFactorisationStrategy()

Solve an affine problem with a single back-substitution per step, refactorising the
tangent only when its [`FactorisationKey`](@ref) changes.
"""
struct CachedFactorisationStrategy <: AbstractSolverStrategy end
