struct StorageResponse{
    Q<:AbstractStorageQuantity,
    D,
    M<:AbstractMeasure,
    L,
} <: AbstractResponse
    quantity::Q
    driving::D
    measure::M
    dependencies::Tuple{Vararg{Symbol}}
    law::L
end

function StorageResponse(;
    quantity::AbstractStorageQuantity,
    rate,
    measure::AbstractMeasure,
    dependencies::Tuple{Vararg{Symbol}}=(),
    law=PlaceholderLaw(:unspecified_storage_law),
)
    return StorageResponse{typeof(quantity),typeof(rate),typeof(measure),typeof(law)}(
        quantity,
        rate,
        measure,
        dependencies,
        law,
    )
end

storage_quantity(response::StorageResponse) = response.quantity
rate_quantity(response::StorageResponse) = response.driving

response_issues(::StorageResponse) = ResponseIssue[]
