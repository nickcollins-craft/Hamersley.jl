#=
Traversal of the dataflow graph.

These helpers walk a vertex and everything beneath it. They are what the formulation verification, the trace
derivation, the linearity inference and the coupling analysis are built on.
=#


"""
    collect_vertices(vertex)

Return a vector holding `vertex` and every vertex beneath it, in depth-first order.
"""
function collect_vertices(vertex::AbstractVertex)
    # We start from an empty vector and let the in-place version append to it recursively.
    vertices = AbstractVertex[]
    collect_vertices!(vertices, vertex)
    return vertices
end


"""
    collect_vertices!(vertices, vertex)

Append `vertex` and every vertex beneath it to `vertices`, in depth-first order.
"""
function collect_vertices!(vertices::Vector{AbstractVertex}, vertex::AbstractVertex)
    # The vertex itself comes before its children, which gives depth-first (pre-)order.
    push!(vertices, vertex)
    for child in subvertices(vertex)
        collect_vertices!(vertices, child)
    end
    return nothing
end


"""
    collect_constitutive(vertex)

Return the distinct [`Constitutive`](@ref) bindings used anywhere beneath `vertex`.
"""
function collect_constitutive(vertex::AbstractVertex)
    bindings = Constitutive[]
    for subvertex in collect_vertices(vertex)
        # Several outputs of the same binding should only report it once.
        if subvertex isa ConstitutiveOutput
            push_unique!(bindings, subvertex.constitutive)
        end
    end
    return bindings
end


"""
    referenced_fields(vertex)

Return every field name that `vertex` or any vertex beneath it refers to, whether or not the
value depends on it (this is what the formulation verification checks for undeclared fields).
"""
function referenced_fields(vertex::AbstractVertex)
    names = Symbol[]
    for subvertex in collect_vertices(vertex)
        # Operators and kinematic vertices name their field; transforms may name a motion.
        if subvertex isa AbstractFieldOperator || subvertex isa AbstractKinematicVertex
            push_unique!(names, operand_field(subvertex))
        elseif subvertex isa AbstractTransformVertex && subvertex.motion !== nothing
            push_unique!(names, subvertex.motion)
        end
    end
    return names
end


"""
    push_unique!(collection, item)

Append `item` to `collection` unless it is already there. Membership is tested with `in`
(that is, `isequal`), which for the immutable vertices, traces and bindings of the graph
compares them by content.
"""
function push_unique!(collection::AbstractVector, item)
    # The collections here hold a handful of items, so a linear search is the clearest.
    if !(item in collection)
        push!(collection, item)
    end
    return nothing
end
