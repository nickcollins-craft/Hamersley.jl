using LinearAlgebra, SparseArrays, Tensors, NonlinearSolve

struct StressResponse{
    S<:AbstractStressMeasure,
    D,
    M<:AbstractMeasure,
    L,
} <: AbstractResponse
    quantity::S
    driving::D
    measure::M
    dependencies::Tuple{Vararg{Symbol}}
    law::L
end

function StressResponse(;
    stress::AbstractStressMeasure,
    conjugate_to,
    measure::AbstractMeasure,
    dependencies::Tuple{Vararg{Symbol}}=(),
    law=PlaceholderLaw(:unspecified_stress_law),
)
    return StressResponse{typeof(stress),typeof(conjugate_to),typeof(measure),typeof(law)}(
        stress,
        conjugate_to,
        measure,
        dependencies,
        law,
    )
end

stress_measure(response::StressResponse) = response.quantity
conjugate_quantity(response::StressResponse) = response.driving

function response_issues(response::StressResponse)
    issues = ResponseIssue[]
    stress = stress_measure(response)
    conjugate_to = conjugate_quantity(response)
    if !is_valid_stress_pair(stress, conjugate_to)
        push!(
            issues,
            ResponseIssue(
                "stress measure $(typeof(stress)) is not work-conjugate with $(typeof(conjugate_to))",
            ),
        )
    end
    if !same_configuration(response_measure(response), conjugate_configuration(conjugate_to))
        push!(
            issues,
            ResponseIssue(
                "stress response measure $(typeof(response_measure(response))) is incompatible with conjugate quantity $(typeof(conjugate_to))",
            ),
        )
    end
    if !same_configuration(response_measure(response), quantity_configuration(stress))
        push!(
            issues,
            ResponseIssue(
                "stress measure $(typeof(stress)) is not expressed over $(typeof(response_measure(response)))",
            ),
        )
    end
    return issues
end
