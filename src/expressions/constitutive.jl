#=
Constitutive nodes.

A `Constitutive` binding connects a law (which declares the tags of its inputs and
outputs) to the nodes that feed each input. Each output of the binding is then a node of
its own, `output(binding, :name)`, which can appear in weak terms, be pulled back, or feed
another law. Separating the law from its binding keeps materials independent of any
particular formulation.
=#


"""
    Constitutive(law, inputs)

The binding of `law` to the nodes that provide its inputs.

# Fields
- `law::L`: the constitutive law.
- `inputs::I`: a named tuple with one node per input slot of `law`.

# Examples
```julia
mechanics = Constitutive(NeoHookean(shear_modulus = 1.0, first_lame_parameter = 2.0),
                         (F = DeformationGradient(:u),))
P = output(mechanics, :P)
```
"""
struct Constitutive{L<:AbstractLaw,I<:NamedTuple}
    law::L
    inputs::I
end


"""
    ConstitutiveOutput(constitutive, name)

The node carrying the output `name` of a [`Constitutive`](@ref) binding. Construct it with
[`output`](@ref).

# Fields
- `constitutive::C`: the binding.
- `name::Symbol`: the name of the output.
"""
struct ConstitutiveOutput{C<:Constitutive} <: AbstractNode
    constitutive::C
    name::Symbol
end


"""
    output(constitutive, name)

Return the node carrying the output `name` of `constitutive`.
"""
output(constitutive::Constitutive, name::Symbol) = ConstitutiveOutput(constitutive, name)


constitutive_law(node::ConstitutiveOutput) = node.constitutive.law


function quantity(node::ConstitutiveOutput, ::Tuple)
    # The output's tag is whatever the law declares for it.
    outputs = law_outputs(constitutive_law(node))
    if haskey(outputs, node.name)
        return outputs[node.name]
    else
        return nothing
    end
end


"""
    dependency_inputs(node)

Return the input nodes that the output `node` depends on, according to the law's
per-output dependencies. Unknown outputs or input slots contribute nothing (the audit
reports them separately).
"""
function dependency_inputs(node::ConstitutiveOutput)
    # We look up which input slots this particular output depends on.
    law = constitutive_law(node)
    dependencies = output_dependencies(law)
    if !haskey(dependencies, node.name)
        return AbstractNode[]
    end

    # We keep only the slots that the binding actually supplies.
    bound_inputs = node.constitutive.inputs
    inputs = AbstractNode[]
    for slot in dependencies[node.name]
        if haskey(bound_inputs, slot)
            push!(inputs, bound_inputs[slot])
        end
    end
    return inputs
end


function node_fields(node::ConstitutiveOutput, assumption::AbstractKinematicAssumption)
    # An output depends on the fields of exactly those inputs it is declared to depend on,
    # which is what keeps structurally zero tangent blocks out of the sparsity pattern.
    fields = Symbol[]
    for input in dependency_inputs(node)
        for field in node_fields(input, assumption)
            if !(field in fields)
                push!(fields, field)
            end
        end
    end
    return fields
end


function is_affine(node::ConstitutiveOutput, assumption::AbstractKinematicAssumption)
    # An output is affine when the law is affine and every input it depends on is affine.
    law_is_affine = is_affine(constitutive_law(node))
    inputs_are_affine = all(input -> is_affine(input, assumption), dependency_inputs(node))
    return law_is_affine && inputs_are_affine
end


subnodes(node::ConstitutiveOutput) = Tuple(values(node.constitutive.inputs))


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

    # Every declared slot must be bound, to a node carrying the expected tag.
    for slot in keys(expected)
        if !haskey(bound, slot)
            push!(issues, "law $(law_name(law)) input $slot is not bound to any node")
            continue
        end
        supplied = quantity(bound[slot], fields)
        if supplied !== nothing && supplied != expected[slot]
            push!(
                issues,
                "law $(law_name(law)) input $slot expects $(expected[slot]) " *
                "but is bound to a node carrying $supplied",
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
