#=
Traversal of the dataflow graph.

These helpers walk a node and everything beneath it. They are what the audit, the trace
derivation, the linearity inference and the coupling analysis are built on.
=#


"""
    collect_nodes(node)

Return a vector holding `node` and every node beneath it, in depth-first order.
"""
function collect_nodes(node::AbstractNode)
    nodes = AbstractNode[]
    collect_nodes!(nodes, node)
    return nodes
end


"""
    collect_nodes!(nodes, node)

Append `node` and every node beneath it to `nodes`, in depth-first order.
"""
function collect_nodes!(nodes::Vector{AbstractNode}, node::AbstractNode)
    push!(nodes, node)
    for child in subnodes(node)
        collect_nodes!(nodes, child)
    end
    return nothing
end


"""
    collect_constitutive(node)

Return the distinct [`Constitutive`](@ref) bindings used anywhere beneath `node`.
"""
function collect_constitutive(node::AbstractNode)
    bindings = Constitutive[]
    for subnode in collect_nodes(node)
        if subnode isa ConstitutiveOutput
            binding = subnode.constitutive
            # Several outputs of the same binding should only report it once.
            if !any(existing -> existing === binding, bindings)
                push!(bindings, binding)
            end
        end
    end
    return bindings
end


"""
    referenced_fields(node)

Return every field name that `node` or any node beneath it refers to, whether or not the
value depends on it (this is what the audit checks for undeclared fields).
"""
function referenced_fields(node::AbstractNode)
    names = Symbol[]
    for subnode in collect_nodes(node)
        # Operators and kinematic nodes name their field; transforms may name a motion.
        if subnode isa AbstractFieldOperator || subnode isa AbstractKinematicNode
            push_unique!(names, operand_field(subnode))
        elseif subnode isa AbstractTransformNode && subnode.motion !== nothing
            push_unique!(names, subnode.motion)
        end
    end
    return names
end


"""
    push_unique!(collection, item)

Append `item` to `collection` unless it is already there.
"""
function push_unique!(collection::AbstractVector, item)
    if !(item in collection)
        push!(collection, item)
    end
    return nothing
end
