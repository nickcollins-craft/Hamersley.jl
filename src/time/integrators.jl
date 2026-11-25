#=
Time integrators.

A `Rate` node in a formulation is an ordinary input to the laws that consume it. What
relates the rate to the unknowns of the current step is the time integrator: nearly every
implicit scheme writes ``\dot{x}_{n+1} = a \, x_{n+1} + (\text{history})``, and the factor
``a = ∂\dot{x}/∂x_{n+1}`` is all the linearisation needs from it. Keeping that factor here,
rather than in the laws, means that the same law works with any integrator, in a steady
solve, or in an explicit scheme.
=#


"""
    AbstractTimeIntegrator

Supertype of time integrators.
"""
abstract type AbstractTimeIntegrator end


"""
    QuasiStatic()

No time integration: rates vanish, and so does their contribution to the tangent.
"""
struct QuasiStatic <: AbstractTimeIntegrator end


"""
    BackwardEuler()

The backward Euler scheme, ``\\dot{x}_{n+1} = (x_{n+1} - x_n)/Δt``.
"""
struct BackwardEuler <: AbstractTimeIntegrator end


"""
    Newmark(; beta, gamma)

The Newmark scheme for second-order problems.

# Fields
- `β::T`: the Newmark parameter ``β`` (``1/4`` for the average-acceleration rule).
- `γ::T`: the Newmark parameter ``γ`` (``1/2`` for the average-acceleration rule).

# References
Newmark, N. M. (1959), A method of computation for structural dynamics, *Journal of the
Engineering Mechanics Division* 85(3), 67–94.
"""
struct Newmark{T<:Real} <: AbstractTimeIntegrator
    β::T
    γ::T
end


function Newmark(; beta::Real, gamma::Real)
    # We promote the two parameters to a common type so that the struct is concrete.
    β, γ = promote(beta, gamma)
    integrator = Newmark(β, γ)
    return integrator
end


"""
    rate_factor(integrator, order, Δt)

Return the factor ``a = ∂x^{(k)}_{n+1} / ∂x_{n+1}`` relating the time derivative of order
`order` (``k = 1`` or ``2``) to the current value under `integrator` with time step `Δt`.

| Integrator       | ``k = 1``          | ``k = 2``             |
|------------------|--------------------|-----------------------|
| `QuasiStatic`    | ``0``              | ``0``                 |
| `BackwardEuler`  | ``1/Δt``           | not defined           |
| `Newmark(β, γ)`  | ``γ/(β Δt)``       | ``1/(β Δt^2)``        |
"""
function rate_factor(::QuasiStatic, order::Integer, Δt::Real)
    # Rates vanish in a quasi-static solve, whatever their order.
    a = zero(Δt)
    return a
end

function rate_factor(::BackwardEuler, order::Integer, Δt::Real)
    if order == 1
        a = 1/Δt
        return a
    else
        throw(
            ArgumentError(
                "BackwardEuler only defines first-order rates, please choose from: " *
                "Newmark(beta = …, gamma = …) for second-order rates",
            ),
        )
    end
end

function rate_factor(integrator::Newmark, order::Integer, Δt::Real)
    β = integrator.β
    γ = integrator.γ
    if order == 1
        # ∂u̇/∂u = γ/(βΔt)
        a = γ/(β*Δt)
        return a
    elseif order == 2
        # ∂ü/∂u = 1/(βΔt²)
        a = 1/(β*Δt^2)
        return a
    else
        throw(ArgumentError("Newmark defines rates of order 1 and 2, got $order"))
    end
end
