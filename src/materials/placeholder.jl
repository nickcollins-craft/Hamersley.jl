#=
Placeholder constitutive laws.

While the numerics are not yet written, formulations still need laws with declared inputs,
outputs and dependencies, so that the formulation verification, the trace derivation and the coupling
analysis can be exercised. `PlaceholderLaw` provides exactly that declaration and nothing
more.
=#


"""
    PlaceholderLaw(name; inputs, outputs, dependencies, affine = false, symmetric_tangent = false)

A constitutive law that declares its interface but cannot be evaluated.

# Fields
- `name::Symbol`: a name for the law, used in messages.
- `inputs::NamedTuple`: input slot names mapped to the quantity tags they expect.
- `outputs::NamedTuple`: output names mapped to the quantity tags they produce.
- `dependencies::NamedTuple`: output names mapped to tuples of the input slots they
  depend on. When omitted, every output depends on every input.
- `affine::Bool`: whether the law is affine in its inputs.
- `symmetric_tangent::Bool`: whether the law's tangent is symmetric (see
  [`tangent_is_symmetric`](@ref)).
"""
struct PlaceholderLaw{I<:NamedTuple,O<:NamedTuple,D<:NamedTuple} <: AbstractLaw
    name::Symbol
    inputs::I
    outputs::O
    dependencies::D
    affine::Bool
    symmetric_tangent::Bool
end


function PlaceholderLaw(
    name::Symbol;
    inputs::NamedTuple,
    outputs::NamedTuple,
    dependencies::Union{NamedTuple,Nothing}=nothing,
    affine::Bool=false,
    symmetric_tangent::Bool=false,
)
    # When no dependencies are given, every output depends on every input.
    if dependencies === nothing
        output_names = keys(outputs)
        all_inputs = map(_ -> keys(inputs), output_names)
        declared_dependencies = NamedTuple{output_names}(all_inputs)
    else
        declared_dependencies = dependencies
    end

    law = PlaceholderLaw(name, inputs, outputs, declared_dependencies, affine, symmetric_tangent)
    return law
end


law_name(law::PlaceholderLaw) = law.name
law_inputs(law::PlaceholderLaw) = law.inputs
law_outputs(law::PlaceholderLaw) = law.outputs
output_dependencies(law::PlaceholderLaw) = law.dependencies
is_affine(law::PlaceholderLaw) = law.affine
tangent_is_symmetric(law::PlaceholderLaw) = law.symmetric_tangent
