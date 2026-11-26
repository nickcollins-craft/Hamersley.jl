#=
Constitutive vertices.

A `Constitutive` binding connects a law (which declares the tags of its inputs and
outputs) to the vertices that feed each input. Each output of the binding is then a vertex of
its own, `output(binding, :name)`, which can appear in weak terms, be pulled back, or feed
another law. Separating the law from its binding keeps materials independent of any
particular formulation.
=#


"""
    Constitutive(law, inputs)

The binding of `law` to the vertices that provide its inputs.

# Fields
- `law::L`: the constitutive law.
- `inputs::I`: a named tuple with one vertex per input slot of `law`.

# Examples
```julia
mechanics = Constitutive(
    NeoHookean(shear_modulus = 1.0, first_lame_parameter = 2.0),
    (F = DeformationGradient(:u),),
)
P = output(mechanics, :P)
```
"""
struct Constitutive{L<:AbstractLaw,I<:NamedTuple}
    law::L
    inputs::I
end


"""
    ConstitutiveOutput(constitutive, name)

The vertex carrying the output `name` of a [`Constitutive`](@ref) binding. Construct it with
[`output`](@ref).

# Fields
- `constitutive::C`: the binding.
- `name::Symbol`: the name of the output.
"""
struct ConstitutiveOutput{C<:Constitutive} <: AbstractVertex
    constitutive::C
    name::Symbol
end


"""
    output(constitutive, name)

Return the vertex carrying the output `name` of `constitutive`.
"""
output(constitutive::Constitutive, name::Symbol) = ConstitutiveOutput(constitutive, name)


constitutive_law(vertex::ConstitutiveOutput) = vertex.constitutive.law


function quantity(vertex::ConstitutiveOutput, ::Tuple)
    # The output's tag is whatever the law declares for it.
    outputs = law_outputs(constitutive_law(vertex))
    if haskey(outputs, vertex.name)
        return outputs[vertex.name]
    else
        return nothing
    end
end


"""
    dependency_inputs(vertex)

Return the input vertices that the output `vertex` depends on, according to the law's
per-output dependencies. Unknown outputs or input slots contribute nothing (the formulation verification
reports them separately).
"""
function dependency_inputs(vertex::ConstitutiveOutput)
    # We look up which input slots this particular output depends on.
    law = constitutive_law(vertex)
    dependencies = output_dependencies(law)
    if !haskey(dependencies, vertex.name)
        return AbstractVertex[]
    end

    # We keep only the slots that the binding actually supplies.
    bound_inputs = vertex.constitutive.inputs
    inputs = AbstractVertex[]
    for slot in dependencies[vertex.name]
        if haskey(bound_inputs, slot)
            push!(inputs, bound_inputs[slot])
        end
    end
    return inputs
end


function vertex_fields(vertex::ConstitutiveOutput, assumption::AbstractKinematicAssumption)
    # An output depends on the fields of exactly those inputs it is declared to depend on,
    # which is what keeps structurally zero tangent blocks out of the sparsity pattern.
    fields = Symbol[]
    for input in dependency_inputs(vertex)
        for field in vertex_fields(input, assumption)
            if !(field in fields)
                push!(fields, field)
            end
        end
    end
    return fields
end


function is_affine(vertex::ConstitutiveOutput, assumption::AbstractKinematicAssumption)
    # An output is affine when the law is affine and every input it depends on is affine.
    law_is_affine = is_affine(constitutive_law(vertex))
    inputs_are_affine =
        all(input -> is_affine(input, assumption), dependency_inputs(vertex))
    return law_is_affine && inputs_are_affine
end


subvertices(vertex::ConstitutiveOutput) = Tuple(values(vertex.constitutive.inputs))


"""
    binding_issues(constitutive, fields)

Return a vector of messages describing every way in which the inputs bound in
`constitutive` fail to match the law's declared input slots.
"""
function binding_issues(constitutive::Constitutive, fields::Tuple)
    issues = String[]
    law = constitutive.law
    expected = law_inputs(law)
    bound = constitutive.inputs

    # Every declared slot must be bound, to a vertex carrying the expected tag.
    for slot in keys(expected)
        if !haskey(bound, slot)
            push!(issues, "law $(law_name(law)) input $slot is not bound to any vertex")
            continue
        end
        supplied = quantity(bound[slot], fields)
        if supplied !== nothing && supplied != expected[slot]
            push!(
                issues,
                "law $(law_name(law)) input $slot expects $(expected[slot]) " *
                    "but is bound to a vertex carrying $supplied",
            )
        end
    end

    # Bindings must not name slots the law does not have.
    for slot in keys(bound)
        if !haskey(expected, slot)
            valid = join(string.(keys(expected)), ", ")
            push!(
                issues,
                "law $(law_name(law)) has no input $slot, please choose from: $valid",
            )
        end
    end
    return issues
end
