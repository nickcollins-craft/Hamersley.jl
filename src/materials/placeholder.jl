#=
Placeholder constitutive laws.

While the numerics are not yet written, formulations still need laws with declared inputs,
outputs and dependencies, so that the formulation verification, the trace derivation and the coupling
analysis can be exercised. `PlaceholderLaw` provides exactly that declaration and nothing
more.
=#


"""
    PlaceholderLaw(name; inputs, outputs, dependencies, affine = false,
                   symmetric_tangent = false, initial_state = nothing)

A constitutive law that declares its interface but cannot be evaluated.

# Fields
- `name::Symbol`: a name for the law, used in messages.
- `inputs::NamedTuple`: input slot names mapped to the quantity tags they expect.
- `outputs::NamedTuple`: output names mapped to the quantity tags they produce.
- `dependencies::NamedTuple`: output names mapped to tuples of the input slots they
  depend on. When omitted, every output depends on every input, and so does any output
  the named tuple leaves out (see [`declared_dependencies`](@ref)).
- `affine::Bool`: whether the law is affine in its inputs.
- `symmetric_tangent::Bool`: whether the law's tangent is symmetric (see
  [`tangent_is_symmetric`](@ref)).
- `state::S`: the initial internal state at a quadrature point, or `nothing` for a law
  without internal variables (see [`initial_state`](@ref)). This is what lets a stateful
  law, such as Maxwell viscoelasticity with ``σ_n`` as its state, be declared before it
  can be evaluated.
"""
struct PlaceholderLaw{I<:NamedTuple,O<:NamedTuple,D<:NamedTuple,S} <: AbstractLaw
    name::Symbol
    inputs::I
    outputs::O
    dependencies::D
    affine::Bool
    symmetric_tangent::Bool
    state::S
end


function PlaceholderLaw(
    name::Symbol;
    inputs::NamedTuple,
    outputs::NamedTuple,
    dependencies::Union{NamedTuple,Nothing}=nothing,
    affine::Bool=false,
    symmetric_tangent::Bool=false,
    initial_state=nothing,
)
    # When no dependencies are given, every output depends on every input.
    if dependencies === nothing
        output_names = keys(outputs)
        all_inputs = map(_ -> keys(inputs), output_names)
        dependency_declaration = NamedTuple{output_names}(all_inputs)
    else
        dependency_declaration = dependencies
    end

    law = PlaceholderLaw(
        name,
        inputs,
        outputs,
        dependency_declaration,
        affine,
        symmetric_tangent,
        initial_state,
    )
    return law
end


law_name(law::PlaceholderLaw) = law.name
law_inputs(law::PlaceholderLaw) = law.inputs
law_outputs(law::PlaceholderLaw) = law.outputs
output_dependencies(law::PlaceholderLaw) = law.dependencies
is_affine(law::PlaceholderLaw) = law.affine
tangent_is_symmetric(law::PlaceholderLaw) = law.symmetric_tangent
initial_state(law::PlaceholderLaw) = law.state
