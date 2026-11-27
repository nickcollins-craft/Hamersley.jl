#=
Boundary conditions as boundary terms.

A natural condition is a boundary integral of the weak form, ``∫ δu \\cdot t \\,
\\mathrm{d}A``, and a Robin condition is one whose flux depends on the solution. Compiling
both into ordinary weak terms means that every analysis of the weak form (linearity,
coupling, the dependence of the measure on the motion, and later assembly) treats them
exactly as it treats a term the user wrote, without special cases.

The measure is what distinguishes a dead load from a follower load. A load measured per
unit reference area is integrated over the reference boundary, which is fixed; a load
measured per unit current area is integrated over the current boundary, whose pull-back
by Nanson's formula depends on the motion, so that under finite strains the term is
neither affine nor independent of the motion.

Conditions are compiled on the discretised fields (after the template's translation), and
only where a facet integral exists: an edge force of a second-gradient continuum would need
an edge measure, which does not yet exist, and is a dead load in any case.
=#


"""
    condition_terms(problem, conditions)

Return the boundary terms that the natural and Robin conditions among `conditions` become,
each paired with the equation tested by its field, as a vector of `(equation, term)`
tuples. Conditions on fields without an equation, or on regions without a facet integral,
give no term (the formulation verification reports the former).
"""
function condition_terms(problem::Problem, conditions)
    pairs = Tuple{Equation,AbstractWeakTerm}[]
    for bc in conditions
        term = condition_term(problem, bc)
        if term === nothing
            continue
        end

        # The term belongs to the equation tested by its (translated) field.
        field = operand_field(term_test(term).vertex)
        equation_position = findfirst(eq -> test_field(eq) === field, problem.equations)
        if equation_position !== nothing
            push!(pairs, (problem.equations[equation_position], term))
        end
    end
    return pairs
end


"""
    condition_term(problem, bc)

Return the boundary term that `bc` becomes, or `nothing` when it does not become one.
"""
condition_term(::Problem, ::AbstractBoundaryCondition) = nothing

function condition_term(problem::Problem, bc::NaturalBC)
    # Only facet regions have a boundary measure to integrate over.
    domain = facet_set_name(bc_region(bc))
    if domain === nothing
        return nothing
    end

    # The load does work against the translated trace of the translated field.
    field, trace = translate_trace(problem.template, bc_field(bc), bc_trace(bc))
    test = Variation(trace_operator(field, trace))
    load = BoundaryLoad(bc_quantity(bc), bc_value(bc), bc_components(bc))
    measure = boundary_measure(quantity_configuration(bc_quantity(bc)))
    term = Term(test, load, measure; on=domain)
    return term
end

function condition_term(problem::Problem, bc::RobinBC)
    # A Robin condition is a facet integral of the value trace.
    domain = facet_set_name(bc_region(bc))
    if domain === nothing
        return nothing
    end

    # Its flux is per unit area of the configuration its measure counts.
    field, _ = translate_trace(problem.template, bc_field(bc), ValueTrace())
    flux = RobinFlux(
        field,
        bc.coefficient,
        bc.reference,
        bc_components(bc),
        configuration(bc.measure),
    )
    term = Term(Variation(Value(field)), flux, bc.measure; on=domain)
    return term
end


"""
    trace_operator(field, trace)

Return the vertex whose variation is the test side of a boundary term against `trace` of
`field`: the value itself for the value trace, and a [`BoundaryTrace`](@ref) otherwise.
"""
trace_operator(field::Symbol, ::ValueTrace) = Value(field)
trace_operator(field::Symbol, trace::AbstractTrace) = BoundaryTrace(field, trace)


"""
    boundary_measure(configuration)

Return the boundary measure that a quantity measured per unit area of `configuration` is
integrated over: the current boundary for a follower load, and the reference boundary
otherwise (including quantities that belong to no configuration).
"""
boundary_measure(::AbstractConfiguration) = ReferenceBoundary()
boundary_measure(::CurrentConfiguration) = CurrentBoundary()
