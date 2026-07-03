struct FluxResponse{
    Q<:AbstractFluxQuantity,
    D,
    M<:AbstractMeasure,
    L,
} <: AbstractResponse
    quantity::Q
    driving::D
    measure::M
    dependencies::Tuple{Vararg{Symbol}}
    law::L
    pulled_back::Bool
end

function FluxResponse(;
    flux::AbstractFluxQuantity,
    driven_by,
    measure::AbstractMeasure,
    dependencies::Tuple{Vararg{Symbol}}=(),
    law=PlaceholderLaw(:unspecified_flux_law),
    pulled_back::Bool=false,
)
    return FluxResponse{typeof(flux),typeof(driven_by),typeof(measure),typeof(law)}(
        flux,
        driven_by,
        measure,
        dependencies,
        law,
        pulled_back,
    )
end

flux_quantity(response::FluxResponse) = response.quantity
driving_operator(response::FluxResponse) = response.driving

function response_issues(response::FluxResponse)
    issues = ResponseIssue[]
    flux = flux_quantity(response)
    driven_by = driving_operator(response)
    flux_config = quantity_configuration(flux)
    driving_config = conjugate_configuration(driven_by)
    measure_config = configuration(response_measure(response))
    if !same_configuration(flux_config, driving_config)
        push!(
            issues,
            ResponseIssue(
                "flux quantity $(typeof(flux)) is incompatible with driving quantity $(typeof(driven_by))",
            ),
        )
    end
    if !same_configuration(flux_config, measure_config) && !response.pulled_back
        push!(
            issues,
            ResponseIssue(
                "flux quantity $(typeof(flux)) cannot be integrated over $(typeof(response_measure(response))) without an explicit pull-back",
            ),
        )
    end
    return issues
end
