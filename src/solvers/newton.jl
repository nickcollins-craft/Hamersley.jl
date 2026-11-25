#=
Newton's method, and the choice of solver strategy.

There is one residual-and-tangent path for every problem. A linear problem is not a
separate assembly path: it is one Newton step from the predictor with a cached
factorisation, which is exact because the residual is affine. What a solver strategy
decides is which tangent blocks to use (all of them for full Newton, only the diagonal
ones for a staggered multiphysics scheme) and whether to reuse a factorisation.
=#


"""
    AbstractSolverStrategy

Supertype of global solution strategies.
"""
abstract type AbstractSolverStrategy end


"""
    AbstractBlockSelection

Supertype of choices of which tangent blocks a strategy assembles and uses. What the graph
can differentiate is kept separate from which blocks a solver actually uses.
"""
abstract type AbstractBlockSelection end


"""
    FullTangent()

Use every block of the tangent: full Newton.
"""
struct FullTangent <: AbstractBlockSelection end


"""
    BlockDiagonalTangent()

Use only the diagonal blocks (each equation differentiated with respect to its own test
field), as a staggered multiphysics scheme does.
"""
struct BlockDiagonalTangent <: AbstractBlockSelection end


"""
    NewtonStrategy(; relative_tolerance, absolute_tolerance, maximum_iterations,
                     block_selection = FullTangent())

Newton's method on the global residual.

# Fields
- `relative_tolerance::T`: convergence tolerance on the residual norm relative to its
  initial value.
- `absolute_tolerance::T`: convergence tolerance on the residual norm itself.
- `maximum_iterations::Int`: the iteration limit before the step is cut back.
- `block_selection::B`: which tangent blocks to use.
"""
struct NewtonStrategy{T<:Real,B<:AbstractBlockSelection} <: AbstractSolverStrategy
    relative_tolerance::T
    absolute_tolerance::T
    maximum_iterations::Int
    block_selection::B
end


function NewtonStrategy(;
    relative_tolerance::Real,
    absolute_tolerance::Real,
    maximum_iterations::Integer,
    block_selection::AbstractBlockSelection = FullTangent(),
)
    # We promote the tolerances to a common type so that the struct is concrete.
    relative, absolute = promote(relative_tolerance, absolute_tolerance)
    strategy = NewtonStrategy(relative, absolute, Int(maximum_iterations), block_selection)
    return strategy
end


"""
    solve_step!(solution, strategy, compiled, state, Δt)

Advance `solution` by one step of size `Δt` with `strategy`, overwriting `solution` and
updating `state`.
"""
function solve_step!(
    solution,
    strategy::AbstractSolverStrategy,
    compiled::CompiledProblem,
    state,
    Δt,
)
    not_yet_implemented("$(nameof(typeof(strategy))) steps")
    return nothing
end


"""
    default_strategy(problem; relative_tolerance, absolute_tolerance, maximum_iterations)

Return the default solver strategy for `problem`: a cached factorisation if its residual
is affine, and full Newton otherwise.
"""
function default_strategy(
    problem::Problem;
    relative_tolerance::Real = 1e-10,
    absolute_tolerance::Real = 1e-12,
    maximum_iterations::Integer = 25,
)
    # The linearity claim has already been checked by the audit, so we can trust it here.
    if is_affine(problem)
        return CachedFactorisationStrategy()
    else
        strategy = NewtonStrategy(;
            relative_tolerance = relative_tolerance,
            absolute_tolerance = absolute_tolerance,
            maximum_iterations = maximum_iterations,
        )
        return strategy
    end
end
