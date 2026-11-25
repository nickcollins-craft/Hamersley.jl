#=
The formulation audit.

Before anything is discretised, we check that a problem makes sense as a formulation:
that its fields exist, that every weak term pairs work-conjugate quantities over the right
configuration, that every law is bound correctly, that boundary conditions target traces
that actually exist, and that a claim of linearity is true. The audit reports every issue
it finds rather than stopping at the first, and it never repairs a problem (for example
by inserting a pull-back): where a fix is obvious, the message suggests it.

The rules, in the order they are applied:

1. Field names are unique, and every field a node refers to is declared.
2. Each term's test side is the variation of a node of the equation's test field.
3. The test side and flux of each term are work-conjugate under the kinematic assumption.
4. The flux of each term is measured per unit of the configuration it is integrated over.
5. Under finite strains, the motion field exists and is a vector field.
6. Each constitutive binding matches its law's declared inputs and outputs.
7. Each field's interpolation family is supported, and element-local fields are reported.
8. Boundary conditions target derived traces, with conjugate quantities, and can be
   imposed; in strict mode, no non-standard trace is left defaulted.
9. A claim that the problem is linear agrees with its structure.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Audit reports ------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    AuditIssue(message)

One problem found by the audit.

# Fields
- `message::String`: a description of the issue, with a suggested fix where one is
  obvious.
"""
struct AuditIssue
    message::String
end


"""
    AuditReport(issues)

The outcome of [`audit`](@ref).

# Fields
- `issues::Vector{AuditIssue}`: every issue found, in rule order.
"""
struct AuditReport
    issues::Vector{AuditIssue}
end


Base.show(io::IO, issue::AuditIssue) = print(io, issue.message)

function Base.show(io::IO, report::AuditReport)
    if isempty(report.issues)
        print(io, "AuditReport: pass")
    else
        print(io, "AuditReport: $(length(report.issues)) issue(s)")
    end
    return nothing
end


"""
    passed(report)

Return `true` when the audit found no issues.
"""
passed(report::AuditReport) = isempty(report.issues)


"""
    issue_messages(report)

Return the messages of every issue in `report`.
"""
issue_messages(report::AuditReport) = String[issue.message for issue in report.issues]


# ---------------------------------------------------------------------------------------- #
# ------------------------------------- The audit ---------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    audit(problem; throw_on_error = false)

Check `problem` against every rule listed at the top of `formulation/audit.jl`.

# Returns
An [`AuditReport`](@ref). With `throw_on_error = true`, an `ArgumentError` listing every
issue is thrown instead when the audit fails.
"""
function audit(problem::Problem; throw_on_error::Bool = false)
    messages = String[]

    # We apply the rules in order, so that the report reads from fields to boundaries.
    append!(messages, audit_fields(problem))
    append!(messages, audit_terms(problem))
    append!(messages, audit_kinematics(problem))
    append!(messages, audit_constitutive(problem))
    append!(messages, audit_interpolations(problem))
    append!(messages, audit_boundary_conditions(problem))
    append!(messages, audit_linearity(problem))

    report = AuditReport(AuditIssue.(messages))
    if throw_on_error && !passed(report)
        throw(ArgumentError(join(messages, "\n")))
    end
    return report
end


# ---------------------------------------------------------------------------------------- #
# -------------------------------- Rule 1: declared fields ------------------------------- #
# ---------------------------------------------------------------------------------------- #

function audit_fields(problem::Problem)
    messages = String[]

    # Field names must be unique.
    seen = Symbol[]
    for name in field_names(problem)
        if name in seen
            push!(messages, "field $name is declared more than once")
        end
        push!(seen, name)
    end

    # Every field that a node refers to must be declared.
    for (equation, term) in problem_terms(problem)
        for node in (term_test(term), term_flux(term))
            for name in referenced_fields(node)
                if !(name in seen)
                    push!(
                        messages,
                        "equation $(equation_name(equation)) refers to undeclared " *
                        "field $name",
                    )
                end
            end
        end
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------- Rules 2–4: test sides, conjugacy, measures ------------------- #
# ---------------------------------------------------------------------------------------- #

function audit_terms(problem::Problem)
    messages = String[]
    assumption = kinematic_assumption(problem)

    for equation in problem.equations
        name = equation_name(equation)

        # Rule 2, first half: the test field itself must be declared.
        if field_by_name(problem, test_field(equation)) === nothing
            push!(
                messages,
                "equation $name uses undeclared test field $(test_field(equation))",
            )
            continue
        end

        for term in equation_terms(equation)
            append!(messages, audit_test_side(problem, equation, term))
            append!(messages, audit_conjugacy(problem, assumption, name, term))
            append!(messages, audit_measure(problem, assumption, name, term))
        end
    end
    return messages
end


"""
    varied_fields(node)

Return the fields that the operators and kinematic nodes beneath `node` act on, ignoring
the motion fields of transforms (which enter the value but not the variation's
direction).
"""
function varied_fields(node::AbstractNode)
    names = Symbol[]
    for subnode in collect_nodes(node)
        if subnode isa AbstractFieldOperator || subnode isa AbstractKinematicNode
            push_unique!(names, operand_field(subnode))
        end
    end
    return names
end


function audit_test_side(problem::Problem, equation::Equation, term::Term)
    messages = String[]
    name = equation_name(equation)
    varied = term_test(term).node

    # Rule 2: the test side must be built from operators, not from constitutive outputs.
    if any(node -> node isa ConstitutiveOutput, collect_nodes(varied))
        push!(
            messages,
            "equation $name varies a constitutive output; the test side of a term must " *
            "be built from field operators and kinematic nodes",
        )
        return messages
    end

    # Rule 2: the varied node must depend on exactly the equation's test field.
    fields = varied_fields(varied)
    if fields != [test_field(equation)]
        push!(
            messages,
            "equation $name is tested by $(test_field(equation)) but one of its terms " *
            "varies $(join(string.(fields), ", "))",
        )
    end
    return messages
end


function audit_conjugacy(
    problem::Problem,
    assumption::AbstractKinematicAssumption,
    name::Symbol,
    term::Term,
)
    messages = String[]
    test_tag = quantity(term_test(term), problem.fields)
    flux_tag = quantity(term_flux(term), problem.fields)

    # A transform whose inner quantity has no counterpart leaves the flux without a tag.
    if flux_tag === nothing && term_flux(term) isa AbstractTransformNode
        inner = quantity(term_flux(term).node, problem.fields)
        push!(messages, "equation $name transforms $inner, which has no counterpart")
        return messages
    end

    # Without both tags (undeclared fields or outputs) there is nothing more to check here.
    if test_tag === nothing || flux_tag === nothing
        return messages
    end

    # Quantities that only make sense under small strains are flagged under finite strains.
    if assumption isa FiniteStrain
        for tag in (test_tag, flux_tag)
            if requires_small_strain(tag)
                push!(
                    messages,
                    "equation $name uses $tag, which is only meaningful under SmallStrain",
                )
            end
        end
    end

    # Rule 3: the two sides of the term must do work on each other.
    if !are_conjugate(assumption, test_tag, flux_tag)
        push!(
            messages,
            "equation $name pairs $test_tag with $flux_tag, which are not work-conjugate " *
            "under $assumption",
        )
    end
    return messages
end


function audit_measure(
    problem::Problem,
    assumption::AbstractKinematicAssumption,
    name::Symbol,
    term::Term,
)
    messages = String[]
    measure = term_measure(term)
    flux_tag = quantity(term_flux(term), problem.fields)

    # Interface terms are part of the structure but cannot yet be assembled.
    if is_interface_term(term)
        push!(
            messages,
            "equation $name has an interface term; integration over interior facets " *
            "is not yet implemented",
        )
    end

    # Rule 4: the flux must be measured per unit of the configuration integrated over.
    if flux_tag !== nothing && !same_configuration(assumption, flux_tag, measure)
        push!(
            messages,
            "equation $name integrates $flux_tag over $measure" *
            measure_hint(flux_tag, measure),
        )
    end
    return messages
end


"""
    measure_hint(flux_tag, measure)

Suggest how to reconcile a flux with the measure it is integrated over. The audit never
inserts transforms itself, because where they go is a modelling decision.
"""
function measure_hint(flux_tag::AbstractQuantity, measure::AbstractMeasure)
    counterpart = pulled_back(flux_tag)
    if configuration(measure) isa ReferenceConfiguration && counterpart !== nothing
        return "; wrap the flux in PullBack(…) to integrate it as $counterpart"
    elseif configuration(measure) isa CurrentConfiguration
        return "; integrate over the matching reference measure, or PushForward(…) the flux"
    else
        return ""
    end
end


# ---------------------------------------------------------------------------------------- #
# --------------------------------- Rule 5: the motion ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

audit_kinematics(::Problem{<:Any,<:Any,<:Any,SmallStrain}) = String[]

function audit_kinematics(problem::Problem{<:Any,<:Any,<:Any,FiniteStrain})
    messages = String[]
    motion = motion_field(kinematic_assumption(problem))
    field = field_by_name(problem, motion)

    # The motion must be a declared vector field, since F = I + ∇_0 u.
    if field === nothing
        push!(messages, "FiniteStrain motion field $motion is not declared")
    elseif !(field_shape(field) isa VectorField)
        push!(messages, "FiniteStrain motion field $motion must be a VectorField")
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------ Rule 6: constitutive bindings --------------------------- #
# ---------------------------------------------------------------------------------------- #

function audit_constitutive(problem::Problem)
    messages = String[]

    # Each binding must match its law's declared input slots.
    for binding in problem_constitutive(problem)
        append!(messages, binding_issues(binding, problem.fields))
        append!(messages, dependency_issues(binding.law))
    end

    # Each output node must name an output the law has.
    for (equation, term) in problem_terms(problem)
        for node in collect_nodes(term_flux(term))
            if node isa ConstitutiveOutput
                outputs = law_outputs(constitutive_law(node))
                if !haskey(outputs, node.name)
                    valid = join(string.(keys(outputs)), ", ")
                    push!(
                        messages,
                        "equation $(equation_name(equation)) uses output $(node.name) of " *
                        "law $(law_name(constitutive_law(node))), please choose from: " *
                        valid,
                    )
                end
            end
        end
    end
    return messages
end


"""
    dependency_issues(law)

Return messages for declared dependencies that name unknown outputs or input slots.
"""
function dependency_issues(law::AbstractLaw)
    messages = String[]
    inputs = law_inputs(law)
    outputs = law_outputs(law)
    for (output_name, slots) in pairs(output_dependencies(law))
        if !haskey(outputs, output_name)
            push!(
                messages,
                "law $(law_name(law)) declares dependencies for unknown output " *
                "$output_name",
            )
        end
        for slot in slots
            if !haskey(inputs, slot)
                push!(
                    messages,
                    "law $(law_name(law)) output $output_name depends on unknown " *
                    "input $slot",
                )
            end
        end
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------ Rule 7: interpolation families -------------------------- #
# ---------------------------------------------------------------------------------------- #

function audit_interpolations(problem::Problem)
    messages = String[]
    for field in problem.fields
        # The single gate on interpolation families.
        spec = field_interpolation(field)
        if !is_supported_family(spec)
            push!(
                messages,
                "field $(field_name(field)) uses $(interpolation_family(spec)), which is " *
                "not yet supported, please choose from: $(supported_family_names())",
            )
        end

        # Static condensation of element-local fields is not yet written.
        if field.element_local
            push!(
                messages,
                "field $(field_name(field)) is element-local; static condensation is " *
                "not yet implemented",
            )
        end
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------ Rule 8: boundary conditions ----------------------------- #
# ---------------------------------------------------------------------------------------- #

function audit_boundary_conditions(problem::Problem)
    messages = String[]
    for bc in problem.boundary_conditions
        append!(messages, audit_boundary_condition(problem, bc))
    end

    # In strict mode, non-standard traces must never be left to the default.
    if problem.strict_boundaries
        for status in boundary_status(problem)
            if status.status isa DefaultedToZero && !is_standard_trace(status.trace)
                push!(
                    messages,
                    "trace $(status.trace) of field $(status.field) on " *
                    "$(status.boundary) " *
                    "is not specified (strict_boundaries = true)",
                )
            end
        end
    end
    return messages
end


function audit_boundary_condition(problem::Problem, bc::AbstractBoundaryCondition)
    messages = String[]

    # The condition's field must be declared.
    if field_by_name(problem, bc_field(bc)) === nothing
        push!(messages, "boundary condition targets undeclared field $(bc_field(bc))")
        return messages
    end

    # We check the condition on the discretised fields, after the template's translation.
    target_name, trace = translate_trace(problem.template, bc_field(bc), bc_trace(bc))
    target = field_by_name(problem, target_name)
    if target === nothing
        push!(
            messages,
            "boundary condition is translated onto undeclared field $target_name",
        )
        return messages
    end

    # The trace must be one that the weak form actually produces for this field.
    derived = derived_traces(problem, target_name)
    if !(trace in derived)
        available = isempty(derived) ? "none" : join(string.(derived), ", ")
        push!(
            messages,
            "trace $trace of field $target_name on $(bc_boundary(bc)) does not appear in " *
            "the weak form (derived traces: $available)",
        )
    end

    # A natural condition's quantity must do work against its trace.
    if bc isa NaturalBC && conjugate_trace(bc_quantity(bc)) != bc_trace(bc)
        push!(
            messages,
            "$(bc_quantity(bc)) is conjugate to $(conjugate_trace(bc_quantity(bc))), " *
            "not to $(bc_trace(bc))",
        )
    end

    # Only supported families have a meaningful conformity to impose against.
    if bc isa EssentialBC && is_supported_family(field_interpolation(target))
        if imposition(bc, target) isa WeakImposition
            push!(
                messages,
                "trace $trace of field $target_name cannot be imposed strongly by its " *
                "interpolation, and weak (Nitsche) imposition is not yet implemented",
            )
        end
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Rule 9: linearity ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

function audit_linearity(problem::Problem)
    messages = String[]

    # The claim is checked once, here, so that the solver can trust it without rechecking.
    if problem.linear === true && !is_affine(problem)
        push!(
            messages,
            "the problem is declared linear but its residual is not affine in the " *
            "unknowns under $(kinematic_assumption(problem))",
        )
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------- Explanation -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    explain(problem)

Return a human-readable summary of `problem`: its kinematic assumption, fields, equations,
laws, whether it is linear, and the status of every derived boundary trace.
"""
function explain(problem::Problem)
    lines = String["Problem under $(kinematic_assumption(problem))"]

    # Fields, with their interpolations.
    push!(lines, "  fields:")
    for field in problem.fields
        spec = field_interpolation(field)
        push!(
            lines,
            "    $(field_name(field)): $(nameof(typeof(field_shape(field)))), " *
            "$(interpolation_family(spec)) order $(interpolation_order(spec))",
        )
    end

    # Equations and laws.
    push!(lines, "  equations: " * join(string.(equation_names(problem)), ", "))
    law_names = [string(law_name(binding.law)) for binding in problem_constitutive(problem)]
    push!(lines, "  laws: " * (isempty(law_names) ? "none" : join(law_names, ", ")))
    push!(lines, "  linear (inferred): $(is_affine(problem))")

    # The status of every derived trace on every named boundary.
    push!(lines, "  boundary traces:")
    for status in boundary_status(problem)
        push!(
            lines,
            "    $(status.boundary): $(status.field) $(status.trace) — " *
            "$(nameof(typeof(status.status)))",
        )
    end
    return join(lines, "\n")
end
