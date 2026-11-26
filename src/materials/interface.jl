#=
The constitutive law interface.

A constitutive law in Hamersley is a map from a named tuple of tagged inputs (kinematic
quantities, field values, their rates) and an old internal state to a named tuple of
tagged outputs and a new internal state. The law declares *what* it consumes and
produces through quantity tags; a `Constitutive` vertex in the formulation (see
`expressions/constitutive.jl`) says *which* vertices feed each input.

Three design decisions are fixed here, ahead of the numerics:

1. A law returns the tangent of its *time-discrete* update, taken with respect to the
   current-step inputs only, with the old state and the time step held fixed. Typically
   this is the algorithmic (consistent) tangent, but a law may return a different one (for
   example a secant or stabilised tangent) when that gives better global behaviour. The
   framework uses whichever tangent the law returns, without checking it against the
   outputs; the outputs alone determine the residual, and so the solution.
2. A law declares, per output, which inputs it depends on. This is what the global
   sparsity pattern and the off-diagonal blocks of a multiphysics tangent are built from.
   The tangent a law returns may be sparser than this declaration, but not denser.
3. A law whose tangent is symmetric says so, so that the corresponding block of the global
   system can be treated as symmetric. This holds for any law deriving from a potential
   ``Ψ`` (whose tangent is a Hessian), but also for a law returning a symmetrised tangent.

Rates of inputs (``\dot{ε}``, ``\dot{T}``) are inputs like any other. The factor that
relates a rate to the current unknowns belongs to the time integrator, never to the law.
=#


"""
    AbstractLaw

Supertype of constitutive laws. A concrete law implements [`law_inputs`](@ref),
[`law_outputs`](@ref) and [`constitutive_update`](@ref), and may refine
[`output_dependencies`](@ref), [`is_affine`](@ref), [`tangent_is_symmetric`](@ref) and
[`initial_state`](@ref).
"""
abstract type AbstractLaw end


"""
    law_name(law)

Return a `Symbol` naming `law`, for messages and reports.
"""
law_name(law::AbstractLaw) = nameof(typeof(law))


"""
    law_inputs(law)

Return a named tuple mapping each input slot of `law` to the quantity tag it expects.
"""
function law_inputs end


"""
    law_outputs(law)

Return a named tuple mapping each output of `law` to the quantity tag it produces.
"""
function law_outputs end


"""
    output_dependencies(law)

Return a named tuple mapping each output of `law` to a tuple of the input slots it
depends on. By default every output depends on every input; a law overrides this to
expose structural zeros in its tangent.
"""
function output_dependencies(law::AbstractLaw)
    # Every output depends on every input unless the law says otherwise.
    input_names = keys(law_inputs(law))
    output_names = keys(law_outputs(law))
    dependencies = NamedTuple{output_names}(map(_ -> input_names, output_names))
    return dependencies
end


"""
    is_affine(law)

Return `true` when every output of `law` is an affine function of its inputs (with
constant coefficients), so that its tangent never changes. Laws are assumed nonlinear
unless they say otherwise.
"""
is_affine(::AbstractLaw) = false


"""
    tangent_is_symmetric(law)

Return `true` when the tangent that `law` returns is symmetric, in the sense that each
block differentiating an output with respect to its work-conjugate input (such as
``∂σ/∂ε`` or ``∂P/∂F``) has major symmetry. This holds for any law deriving from a
potential ``Ψ``, whose tangent is a Hessian, and for any law returning a symmetrised
tangent. Blocks between an output and a non-conjugate input are not covered. Laws are
assumed to have non-symmetric tangents unless they say otherwise.
"""
tangent_is_symmetric(::AbstractLaw) = false


"""
    initial_state(law)

Return the initial internal state of `law` at a quadrature point, or `nothing` for laws
without internal variables.
"""
initial_state(::AbstractLaw) = nothing


"""
    constitutive_update(law, inputs, state_old, Δt)

Evaluate the time-discrete update of `law`.

# Arguments
- `law`: the constitutive law, holding its material parameters.
- `inputs`: a named tuple with one value per input slot of `law`, at the current step.
- `state_old`: the converged internal state at the previous step.
- `Δt`: the time step.

# Returns
A tuple `(outputs, tangents, state_new)`, where `outputs` has one value per output of the
law, `tangents` holds the derivative of each output with respect to each input it depends
on, and `state_new` is the updated internal state. The tangents are typically the
algorithmic (consistent) tangent of the update, but a law may return a different one; the
framework uses them as given.
"""
function constitutive_update(law::AbstractLaw, inputs, state_old, Δt)
    not_yet_implemented("constitutive_update for $(law_name(law))")
    return nothing
end
