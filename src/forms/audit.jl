struct AuditIssue
    message::String
end

struct AuditReport
    issues::Vector{AuditIssue}
end

Base.show(io::IO, issue::AuditIssue) = print(io, issue.message)
Base.show(io::IO, report::AuditReport) =
    isempty(report.issues) ? print(io, "AuditReport: pass") : print(io, "AuditReport: $(length(report.issues)) issue(s)")

passed(report::AuditReport) = isempty(report.issues)

function audit(problem::Problem; throw_on_error::Bool=false)
    issues = AuditIssue[]
    append!(issues, audit_fields(problem))
    append!(issues, audit_equations(problem))
    append!(issues, audit_responses(problem))
    append!(issues, audit_boundary_conditions(problem))
    report = AuditReport(issues)
    if throw_on_error && !passed(report)
        throw(ArgumentError(join((issue.message for issue in issues), "\n")))
    end
    return report
end

function explain(problem::Problem)
    lines = String[
        "Problem",
        "  fields: " * join(string.(field_names(problem)), ", "),
        "  equations: " * join(string.(equation_names(problem)), ", "),
        "  responses: " * string(length(problem_responses(problem))),
        "  boundary conditions: " * string(length(problem_boundary_conditions(problem))),
    ]
    return join(lines, "\n")
end

function audit_fields(problem::Problem)
    issues = AuditIssue[]
    names = Symbol[]
    for field in problem_fields(problem)
        name = field_name(field)
        if name in names
            push!(issues, AuditIssue("field $name is declared more than once"))
        end
        push!(names, name)
    end
    return issues
end

function audit_equations(problem::Problem)
    issues = AuditIssue[]
    for equation in problem_equations(problem)
        if field_by_name(problem, test_field(equation)) === nothing
            push!(issues, AuditIssue("equation $(equation_name(equation)) uses undeclared test field $(test_field(equation))"))
        end
    end
    return issues
end

function audit_responses(problem::Problem)
    issues = AuditIssue[]
    known = Set(field_names(problem))
    for response in problem_responses(problem)
        for issue in response_issues(response)
            push!(issues, AuditIssue(issue.message))
        end
        for dependency in response_dependencies(response)
            if !(dependency in known)
                push!(issues, AuditIssue("response depends on undeclared field $dependency"))
            end
        end
    end
    return issues
end

function audit_boundary_conditions(problem::Problem)
    issues = AuditIssue[]
    for bc in problem_boundary_conditions(problem)
        field = field_by_name(problem, bc_field(bc))
        if field === nothing
            push!(issues, AuditIssue("boundary condition targets undeclared field $(bc_field(bc))"))
            continue
        end
        trace = bc_trace(bc)
        if trace === nothing
            push!(issues, AuditIssue("natural boundary condition on $(bc_field(bc)) must target an explicit trace"))
            continue
        end
        if !supports_trace(problem, bc_field(bc), trace)
            push!(issues, AuditIssue("trace $(typeof(trace)) is not supported for field $(bc_field(bc))"))
        end
        if bc isa EssentialBC && !supports_strong_trace(field_space(field), trace)
            push!(issues, AuditIssue("field $(bc_field(bc)) cannot strongly prescribe trace $(typeof(trace)) with $(typeof(field_space(field)))"))
        end
        if bc isa NaturalBC && !same_configuration(quantity_configuration(bc_quantity(bc)), field_configuration(field))
            push!(issues, AuditIssue("natural boundary quantity $(typeof(bc_quantity(bc))) is incompatible with field $(bc_field(bc)) configuration"))
        end
    end
    return issues
end
