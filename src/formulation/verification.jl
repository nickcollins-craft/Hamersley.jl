#=
The formulation verification.

Before anything is discretised, we check that a problem makes sense as a formulation:
that its fields exist, that every weak term pairs work-conjugate quantities over the right
configuration, that every law is bound correctly, that boundary conditions target traces
that actually exist, and that a claim of linearity is true. The formulation verification reports every issue
it finds rather than stopping at the first, and it never repairs a problem (for example
by inserting a pull-back): where a fix is obvious, the message suggests it.

The rules, in the order they are applied:

1. Field names are unique, and every field a vertex refers to is declared.
2. Each term's test side is the variation of a vertex of the equation's test field.
3. The test side and flux of each term are work-conjugate under the kinematic assumption.
4. The flux of each term is measured per unit of the configuration it is integrated over.
5. Under finite strains, the motion field exists and is a vector field.
6. Each constitutive binding matches its law's declared inputs and outputs.
7. Each field's interpolation family is supported, and element-local fields are reported.
8. Boundary conditions target derived traces, with conjugate quantities, and can be
   imposed; in strict mode, no non-standard trace is left defaulted.
9. A claim that the problem is linear agrees with its structure.
10. Every equation's residual depends on at least one unknown, since otherwise its rows of
    the tangent are empty and the system is singular.
=#


# ---------------------------------------------------------------------------------------- #
# ------------------------------------ Verification reports ------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    VerificationIssue(message)

One problem found by the verification.

# Fields
- `message::String`: a description of the issue, with a suggested fix where one is
  obvious.
"""
struct VerificationIssue
    message::String
end


"""
    VerificationReport(issues)

The outcome of [`verification`](@ref).

# Fields
- `issues::Vector{VerificationIssue}`: every issue found, in rule order.
"""
struct VerificationReport
    issues::Vector{VerificationIssue}
end


Base.show(io::IO, issue::VerificationIssue) = print(io, issue.message)

function Base.show(io::IO, report::VerificationReport)
    # We print a one-line summary; the messages themselves are in `issue_messages`.
    if isempty(report.issues)
        print(io, "VerificationReport: pass")
    else
        print(io, "VerificationReport: $(length(report.issues)) issue(s)")
    end
    return nothing
end


"""
    passed(report)

Return `true` when the verification found no issues.
"""
passed(report::VerificationReport) = isempty(report.issues)


"""
    issue_messages(report)

Return the messages of every issue in `report`.
"""
issue_messages(report::VerificationReport) =
    String[issue.message for issue in report.issues]


# ---------------------------------------------------------------------------------------- #
# -------------------------------------- The verification --------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification(problem; throw_on_error = false)

Check `problem` against every rule listed at the top of `formulation/verification.jl`.

# Returns
An [`VerificationReport`](@ref). With `throw_on_error = true`, an `ArgumentError` listing every
issue is thrown instead when the verification fails.
"""
function verification(problem::Problem; throw_on_error::Bool=false)
    messages = String[]

    # We apply the rules in order, so that the report reads from fields to boundaries.
    append!(messages, verification_fields(problem))
    append!(messages, verification_terms(problem))
    append!(messages, verification_kinematics(problem))
    append!(messages, verification_constitutive(problem))
    append!(messages, verification_interpolations(problem))
    append!(messages, verification_boundary_conditions(problem))
    append!(messages, verification_linearity(problem))
    append!(messages, verification_dependence(problem))

    report = VerificationReport(VerificationIssue.(messages))
    if throw_on_error && !passed(report)
        throw(ArgumentError(join(messages, "\n")))
    end
    return report
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------- Rule 1: declared fields -------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_fields(problem)

Rule 1: return a message for every field name declared more than once, and for every
field that a vertex refers to but `problem` does not declare.
"""
function verification_fields(problem::Problem)
    messages = String[]

    # Field names must be unique.
    seen = Symbol[]
    for name in field_names(problem)
        if name in seen
            push!(messages, "field $name is declared more than once")
        end
        push!(seen, name)
    end

    # Every field that a vertex refers to must be declared.
    for (equation, term) in problem_terms(problem)
        for vertex in (term_test(term), term_flux(term))
            for name in referenced_fields(vertex)
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
# ---------------------- Rules 2–4: test sides, conjugacy, measures ---------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_terms(problem)

Rules 2 to 4: return a message for every equation whose test field is undeclared, and for
every term whose test side, conjugacy or measure is wrong (see [`verification_test_side`](@ref),
[`verification_conjugacy`](@ref) and [`verification_measure`](@ref)).
"""
function verification_terms(problem::Problem)
    messages = String[]
    assumption = kinematic_assumption(problem)

    # We verification equation by equation, so that every message can name its equation.
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

        # Rules 2 (second half) to 4 apply to each term in turn.
        for term in equation_terms(equation)
            append!(messages, verification_test_side(problem, equation, term))
            append!(messages, verification_conjugacy(problem, assumption, name, term))
            append!(messages, verification_measure(problem, assumption, name, term))
        end
    end
    return messages
end


"""
    varied_fields(vertex)

Return the fields that the operators and kinematic vertices beneath `vertex` act on, ignoring
the motion fields of transforms (which enter the value but not the variation's
direction).
"""
function varied_fields(vertex::AbstractVertex)
    names = Symbol[]
    for subvertex in collect_vertices(vertex)
        if subvertex isa AbstractFieldOperator || subvertex isa AbstractKinematicVertex
            push_unique!(names, operand_field(subvertex))
        end
    end
    return names
end


"""
    verification_test_side(problem, equation, term)

Rule 2: return a message if the test side of `term` is not the variation of a vertex built
from field operators and kinematic vertices on the test field of `equation` alone.
"""
function verification_test_side(problem::Problem, equation::Equation, term::Term)
    messages = String[]
    name = equation_name(equation)
    varied = term_test(term).vertex

    # Rule 2: the test side must be built from operators, not from constitutive outputs.
    if any(vertex -> vertex isa ConstitutiveOutput, collect_vertices(varied))
        push!(
            messages,
            "equation $name varies a constitutive output; the test side of a term must " *
                "be built from field operators and kinematic vertices",
        )
        return messages
    end

    # Rule 2: the varied vertex must depend on exactly the equation's test field.
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


"""
    verification_conjugacy(problem, assumption, name, term)

Rule 3: return a message if the test side and flux of `term` (in the equation called
`name`) are not work-conjugate under `assumption`, or if either is a quantity that is only
meaningful under small strains while `assumption` is finite.
"""
function verification_conjugacy(
    problem::Problem,
    assumption::AbstractKinematicAssumption,
    name::Symbol,
    term::Term,
)
    # We first find the quantity tags of the two sides of the term.
    messages = String[]
    test_tag = quantity(term_test(term), problem.fields)
    flux_tag = quantity(term_flux(term), problem.fields)

    # A transform whose inner quantity has no counterpart leaves the flux without a tag.
    if flux_tag === nothing && term_flux(term) isa AbstractTransformVertex
        inner = quantity(term_flux(term).vertex, problem.fields)
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


"""
    verification_measure(problem, assumption, name, term)

Rule 4: return a message if the flux of `term` (in the equation called `name`) is not
measured per unit of the configuration it is integrated over, suggesting the transform
that would reconcile them, and report interface terms, which cannot yet be assembled.
"""
function verification_measure(
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

Suggest how to reconcile a flux with the measure it is integrated over. The formulation verification never
inserts transforms itself, because where they go is a modelling decision.
"""
function measure_hint(flux_tag::AbstractQuantity, measure::AbstractMeasure)
    # Over a reference measure, a spatial flux with a reference counterpart can be pulled
    # back; over a current measure, either the measure or the flux has to change.
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
# ---------------------------------- Rule 5: the motion ---------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_kinematics(problem)

Rule 5: under finite strains, return a message if the motion field is undeclared or is
not a vector field. Under small strains there is no motion to check.
"""
verification_kinematics(::Problem{<:Any,<:Any,<:Any,SmallStrain}) = String[]

function verification_kinematics(problem::Problem{<:Any,<:Any,<:Any,FiniteStrain})
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
# ---------------------------- Rule 6: constitutive bindings ----------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_constitutive(problem)

Rule 6: return a message for every binding whose inputs do not match its law's declared
slots, every inconsistent dependency declaration, and every output vertex that names an
output its law does not have.
"""
function verification_constitutive(problem::Problem)
    messages = String[]

    # Each binding must match its law's declared input slots.
    for binding in problem_constitutive(problem)
        append!(messages, binding_issues(binding, problem.fields))
        append!(messages, dependency_issues(binding.law))
    end

    # Each output vertex must name an output the law has.
    for (equation, term) in problem_terms(problem)
        for vertex in collect_vertices(term_flux(term))
            if vertex isa ConstitutiveOutput
                outputs = law_outputs(constitutive_law(vertex))
                if !haskey(outputs, vertex.name)
                    valid = join(string.(keys(outputs)), ", ")
                    push!(
                        messages,
                        "equation $(equation_name(equation)) uses output $(vertex.name) of " *
                            "law $(law_name(constitutive_law(vertex))), please choose from: " *
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
        # Each declared output must be one that the law actually produces.
        if !haskey(outputs, output_name)
            push!(
                messages,
                "law $(law_name(law)) declares dependencies for unknown output " *
                    "$output_name",
            )
        end

        # Each slot that the output depends on must be one of the law's inputs.
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
# ---------------------------- Rule 7: interpolation families ---------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_interpolations(problem)

Rule 7: return a message for every field whose interpolation family is not supported, and
for every element-local field (whose static condensation is not yet written).
"""
function verification_interpolations(problem::Problem)
    messages = String[]
    for field in problem.fields
        # The single gate on interpolation families.
        if !is_supported_family(field_interpolation(field))
            push!(messages, unsupported_family_message(field))
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
# ----------------------------- Rule 8: boundary conditions ------------------------------ #
# ---------------------------------------------------------------------------------------- #

"""
    verification_boundary_conditions(problem)

Rule 8: return the messages of [`verification_schedule`](@ref), of
[`verification_boundary_condition`](@ref) for every condition as written, and, phase by
phase, of [`verification_conflicts`](@ref) and (in strict mode) one for every non-standard
trace that is left defaulted to zero. When the problem has several phases, the messages
that belong to one phase name it.
"""
function verification_boundary_conditions(problem::Problem)
    messages = String[]

    # The schedule itself must make sense.
    append!(messages, verification_schedule(problem))

    # Each condition is checked once, as it was written.
    for bc in written_conditions(problem)
        append!(messages, verification_boundary_condition(problem, bc))
    end

    # Conflicts and defaulted traces depend on the phase, after inheritance.
    phases = problem_phases(problem)
    for (index, phase) in enumerate(phases)
        prefix = length(phases) > 1 ? "in phase $(phase_name(phase)), " : ""
        conditions = phase_conditions(problem, index)
        for message in verification_conflicts(problem, conditions)
            push!(messages, prefix * message)
        end
        if problem.strict_boundaries
            for message in strict_boundary_issues(problem, index)
                push!(messages, prefix * message)
            end
        end
    end
    return messages
end


"""
    strict_boundary_issues(problem, phase)

Return a message for every non-standard trace that no claim covers at all in `phase`.
Traces covered on only some components leave a remainder whose existence depends on the
spatial dimension, which the compile step checks.
"""
function strict_boundary_issues(problem::Problem, phase)
    messages = String[]
    for status in boundary_status(problem, phase)
        if status.status isa DefaultedToZero &&
           selects_all(status.components) &&
           !is_standard_trace(status.trace)
            push!(
                messages,
                "trace $(status.trace) of field $(status.field) on " *
                    "$(region_label(status.region)) is not specified " *
                    "(strict_boundaries = true)",
            )
        end
    end
    return messages
end


"""
    verification_schedule(problem)

Return a message for every phase name used twice, every phase that can never be reached
(after a phase that never ends, or ending no later than the phase before it), and every
release that removes nothing.
"""
function verification_schedule(problem::Problem)
    messages = String[]
    phases = problem_phases(problem)

    # Phase names must be unique, since messages and reports name phases.
    seen = Symbol[]
    for phase in phases
        if phase_name(phase) in seen
            push!(messages, "phase $(phase_name(phase)) is declared more than once")
        end
        push!(seen, phase_name(phase))
    end

    # Phases must end in order, and only the last may never end.
    previous_time = 0.0
    for (index, phase) in enumerate(phases)
        if phase.until isa Never && index < length(phases)
            push!(
                messages,
                "phase $(phase_name(phase)) never ends, so phase " *
                    "$(phase_name(phases[index+1])) is never reached",
            )
        elseif phase.until isa AtTime
            if phase.until.time <= previous_time
                push!(
                    messages,
                    "phase $(phase_name(phase)) ends at time $(phase.until.time), " *
                        "which is not after $(previous_time)",
                )
            end
            previous_time = phase.until.time
        end
    end

    # A release that removes nothing is almost always a mistyped region or field.
    for (phase, release) in ineffective_releases(problem.schedule)
        push!(
            messages,
            "phase $(phase_name(phase)) releases trace $(bc_trace(release)) of field " *
                "$(bc_field(release)) on $(region_label(bc_region(release)))" *
                "$(selection_label(bc_components(release))), which it does not inherit",
        )
    end
    return messages
end


"""
    verification_boundary_condition(problem, bc)

Return a message if `bc` targets an undeclared field; acts on a kind of region it cannot
act on; targets a trace that the weak form does not produce; pairs a quantity with a trace
it is not conjugate to; cannot be imposed by the field's interpolation; selects components
badly; has a value of the wrong shape; depends on the current position without being a
follower load; or (for a Robin condition) has no flux for its field.
"""
function verification_boundary_condition(problem::Problem, bc::AbstractBoundaryCondition)
    messages = String[]
    location = "on $(region_label(bc_region(bc)))"

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

    # The condition must act on a kind of region that can carry it.
    issue = region_issue(bc)
    if issue !== nothing
        push!(messages, "boundary condition on field $(bc_field(bc)): " * issue)
    end

    # The trace must be one that the weak form actually produces for this field; values at
    # points or along edges are discrete, and are not traces of the weak form.
    derived = derived_traces(problem, target_name)
    if !is_discrete_condition(bc) && !(trace in derived)
        available = isempty(derived) ? "none" : join(string.(derived), ", ")
        push!(
            messages,
            "trace $trace of field $target_name $location does not appear in " *
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

    # The components, the value and its dependencies are checked against the field.
    append!(messages, verification_condition_value(bc, target, location))
    return messages
end


"""
    verification_condition_value(bc, target, location)

Return a message for every badly formed component selection of `bc`, every value of the
wrong shape for the field `target`, a dependence on the current position by anything but
a follower load, and (for a Robin condition) a field without a Robin flux or a measure
that is not a boundary measure. `location` names the region in the messages.
"""
function verification_condition_value(
    bc::AbstractBoundaryCondition,
    target::Field,
    location::String,
)
    messages = String[]
    label = "$(condition_label(bc)) on field $(bc_field(bc)) $location"
    components = bc_components(bc)
    full_shape = full_value_shape(bc, field_shape(target))

    # A selection must be well formed, and a scalar has no components to select.
    for issue in selection_issues(components)
        push!(messages, "$label $issue")
    end
    if full_shape isa ScalarValue && !selects_all(components)
        push!(messages, "$label selects components of a scalar")
    end

    # The value (or the reference value of a Robin condition) must have the right shape.
    shape = selected_shape(full_shape, components)
    for issue in prescription_shape_issues(bc_value(bc), shape)
        push!(messages, "value of $label $issue")
    end

    # Only a follower load may depend on the current position.
    is_follower = bc isa NaturalBC && is_follower_quantity(bc_quantity(bc))
    if depends_on_current_position(bc_value(bc)) && !is_follower
        push!(
            messages,
            "value of $label depends on the current position, which only a follower " *
                "load (a quantity per unit current area) may",
        )
    end

    # A Robin condition needs a flux for its field, and a boundary measure.
    if bc isa RobinBC
        if robin_quantity(field_shape(target), configuration(bc.measure)) === nothing
            push!(messages, "$label has no Robin flux for a field of its shape")
        end
        if !(bc.measure isa AbstractBoundaryMeasure)
            push!(messages, "$label is integrated over $(bc.measure), not a boundary")
        end
    end
    return messages
end


"""
    verification_conflicts(problem, conditions)

Return a message for every pair of claims (by `conditions`, or by boundary terms written in
the equations) on overlapping components of the same trace on the same region that cannot
both hold: anything together with an essential condition, two natural conditions, or two
Robin conditions. A natural condition, a Robin condition and boundary terms may share
components, since they are all boundary integrals and their sum is explicit.
"""
function verification_conflicts(problem::Problem, conditions)
    messages = String[]
    claims = trace_claims(problem, conditions; include_discrete=true)
    for i in eachindex(claims), j in eachindex(claims)
        a = claims[i]
        b = claims[j]

        # Each pair is considered once, and only on the same trace of the same region.
        if j <= i || a.region != b.region || a.field !== b.field || a.trace != b.trace
            continue
        end
        if !overlaps(a.components, b.components) || !claims_conflict(a, b)
            continue
        end
        push!(
            messages,
            "trace $(a.trace) of field $(a.field) on $(region_label(a.region)) is " *
                "claimed by both $(claim_label(a)) and $(claim_label(b))",
        )
    end
    return messages
end


"""
    claims_conflict(a, b)

Return `true` when the claims `a` and `b` on the same components cannot both hold.
"""
function claims_conflict(a::TraceClaim, b::TraceClaim)
    # An essential condition fixes the trace, so nothing else can claim it.
    if a.status isa PrescribedTrace || b.status isa PrescribedTrace
        return true
    end

    # Two conditions of the same natural kind would be silently superposed.
    both_natural = a.source isa NaturalBC && b.source isa NaturalBC
    both_robin = a.source isa RobinBC && b.source isa RobinBC
    return both_natural || both_robin
end


"""
    claim_label(claim)

Return a description of what makes `claim`, for messages.
"""
function claim_label(claim::TraceClaim)
    label = if claim.source isa AbstractBoundaryCondition
        condition_label(claim.source)
    else
        "a boundary term"
    end
    return label * selection_label(claim.components)
end


# ---------------------------------------------------------------------------------------- #
# ---------------------------------- Rule 9: linearity ----------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_linearity(problem)

Rule 9: return a message if `problem` is declared linear but its residual is not affine
in the unknowns.
"""
function verification_linearity(problem::Problem)
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
# ------------------------------ Rule 10: dependence on unknowns ------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    verification_dependence(problem)

Rule 10: return a message for every equation whose residual depends on none of the
unknowns, naming the laws whose dependency declarations should be checked.

Such an equation contributes empty rows to the tangent, so the system is singular whatever
the formulation. Saddle points are unaffected: the incompressibility equation tested by
``δp`` has a zero diagonal block, but still depends on ``u``.
"""
function verification_dependence(problem::Problem)
    messages = String[]
    assumption = kinematic_assumption(problem)
    unknowns = field_names(problem)

    for equation in problem.equations
        name = equation_name(equation)

        # An undeclared test field has already been reported by rule 2, and a flux without
        # a tag (an undeclared field or an unknown output) by rules 1, 3 or 6, so we do
        # not repeat them as a missing dependence.
        if field_by_name(problem, test_field(equation)) === nothing
            continue
        end
        untagged = term -> quantity(term_flux(term), problem.fields) === nothing
        if any(untagged, equation_terms(equation))
            continue
        end

        # An equation without terms has nothing to depend on.
        terms = equation_terms(equation)
        if isempty(terms)
            push!(
                messages,
                "equation $name has no terms, so its rows of the tangent are empty",
            )
            continue
        end

        # We gather every unknown the equation's terms depend on, and the laws they use.
        dependencies = Symbol[]
        laws = Symbol[]
        for term in terms
            for field in term_dependencies(term, assumption)
                push_unique!(dependencies, field)
            end
            for binding in collect_constitutive(term_flux(term))
                push_unique!(laws, law_name(binding.law))
            end
        end

        # Undeclared fields were reported by rule 1, so only declared ones count here.
        if !any(field -> field in unknowns, dependencies)
            hint = if isempty(laws)
                ""
            else
                "; check the dependency declarations of " * join(string.(laws), ", ")
            end
            push!(
                messages,
                "equation $name does not depend on any unknown, so its rows of the " *
                    "tangent are empty" * hint,
            )
        end
    end
    return messages
end


# ---------------------------------------------------------------------------------------- #
# ------------------------------------- Explanation -------------------------------------- #
# ---------------------------------------------------------------------------------------- #

"""
    explain(problem)

Return a human-readable summary of `problem`: its kinematic assumption, fields, equations,
laws, whether it is linear, and, phase by phase, the status of every derived boundary trace
and the discrete conditions (at points and along edges).
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

    # The status of every derived trace on every named region, phase by phase.
    phases = problem_phases(problem)
    for (index, phase) in enumerate(phases)
        indent = "  "
        if length(phases) > 1
            push!(lines, "  phase $(phase_name(phase)) (until $(phase.until)):")
            indent = "    "
        end
        append!(lines, explain_phase(problem, index, indent))
    end
    return join(lines, "\n")
end


"""
    explain_phase(problem, phase, indent)

Return the lines of [`explain`](@ref) for the phase `phase`: the status of every derived
trace, and the discrete conditions, indented by `indent`.
"""
function explain_phase(problem::Problem, phase, indent::String)
    lines = String[]

    # Every derived trace, component by component.
    push!(lines, indent * "boundary traces:")
    for status in boundary_status(problem, phase)
        push!(
            lines,
            indent *
            "  $(region_label(status.region)): $(status.field) $(status.trace)" *
            "$(selection_label(status.components)) — $(status_label(status.status))",
        )
    end

    # Conditions at points and along edges are not traces, so they are listed apart.
    discrete = filter(is_discrete_condition, phase_conditions(problem, phase))
    if !isempty(discrete)
        push!(lines, indent * "discrete conditions:")
        for bc in discrete
            push!(
                lines,
                indent *
                "  $(region_label(bc_region(bc))): $(bc_field(bc)) $(bc_trace(bc))" *
                "$(selection_label(bc_components(bc))) — Prescribed",
            )
        end
    end
    return lines
end
