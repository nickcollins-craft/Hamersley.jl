#=
Field coupling.

Which blocks of the global tangent can be non-zero follows from the graph: the equation
tested by field ``i`` couples to field ``j`` exactly when some term of that equation
depends on ``j``. Because laws declare their dependencies per output, structurally zero
blocks (a heat flux that does not depend on the pore pressure, say) stay out of the
sparsity pattern. The matrix built here is in the form Ferrite's `add_sparsity_entries!`
accepts through its `coupling` keyword: one row and one column per field, in the order the
problem declares them.
=#


"""
    coupling_matrix(problem)

Return a `Matrix{Bool}` whose entry `(i, j)` is `true` when the equation tested by field
`i` depends on field `j`.
"""
function coupling_matrix(problem::Problem)
    names = field_names(problem)
    number_of_fields = length(names)
    coupling = falses(number_of_fields, number_of_fields)
    assumption = kinematic_assumption(problem)

    # Every term couples its test field to every field its flux depends on.
    for (equation, term) in problem_terms(problem)
        row = findfirst(==(test_field(equation)), names)
        if row === nothing
            continue
        end
        mark_dependencies!(coupling, row, vertex_fields(term_flux(term), assumption), names)

        # The test side only adds dependencies when its variation depends on the solution
        # (δE = sym(Fᵀ ∇_0 δu) does; δp and δF = ∇_0 δu do not).
        test = term_test(term)
        if !is_affine(test, assumption)
            mark_dependencies!(coupling, row, vertex_fields(test, assumption), names)
        end
    end

    # Follower loads couple the loaded field to the motion.
    motion = motion_field(assumption)
    for bc in problem.boundary_conditions
        row = findfirst(==(bc_field(bc)), names)
        if row !== nothing && motion !== nothing && is_follower_load(bc, assumption)
            mark_dependencies!(coupling, row, [motion], names)
        end
    end
    return Matrix{Bool}(coupling)
end


"""
    mark_dependencies!(coupling, row, dependencies, names)

Set `coupling[row, j]` to `true` for every field in `dependencies`, where `j` is that
field's position in `names`. Undeclared fields are skipped (the formulation verification reports them).
"""
function mark_dependencies!(
    coupling::AbstractMatrix{Bool},
    row::Integer,
    dependencies::AbstractVector{Symbol},
    names::AbstractVector{Symbol},
)
    # The column of each field is its position in the declaration order.
    for dependency in dependencies
        column = findfirst(==(dependency), names)
        if column !== nothing
            coupling[row, column] = true
        end
    end
    return nothing
end
