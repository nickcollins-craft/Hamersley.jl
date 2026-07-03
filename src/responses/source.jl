struct SourceResponse{
    Q<:AbstractSourceQuantity,
    M<:AbstractMeasure,
    L,
} <: AbstractResponse
    quantity::Q
    driving::Nothing
    measure::M
    dependencies::Tuple{Vararg{Symbol}}
    law::L
end

function SourceResponse(;
    source::AbstractSourceQuantity,
    measure::AbstractMeasure,
    dependencies::Tuple{Vararg{Symbol}}=(),
    law=PlaceholderLaw(:unspecified_source_law),
)
    return SourceResponse{typeof(source),typeof(measure),typeof(law)}(
        source,
        nothing,
        measure,
        dependencies,
        law,
    )
end

source_quantity(response::SourceResponse) = response.quantity

response_issues(::SourceResponse) = ResponseIssue[]
