using Ferrite

struct FerriteBackend end

function lower_to_backend(::FerriteBackend, problem::Problem)
    report = audit(problem)
    passed(report) || throw(ArgumentError("problem failed audit before Ferrite lowering"))
    return problem
end
