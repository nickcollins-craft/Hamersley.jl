abstract type AbstractConfiguration end

struct ReferenceConfiguration <: AbstractConfiguration end
struct CurrentConfiguration <: AbstractConfiguration end
struct IntermediateConfiguration <: AbstractConfiguration end

abstract type AbstractMeasure end
abstract type AbstractVolumeMeasure <: AbstractMeasure end
abstract type AbstractBoundaryMeasure <: AbstractMeasure end

struct ReferenceVolume <: AbstractVolumeMeasure end
struct CurrentVolume <: AbstractVolumeMeasure end
struct ReferenceBoundary <: AbstractBoundaryMeasure end
struct CurrentBoundary <: AbstractBoundaryMeasure end

configuration(::ReferenceVolume) = ReferenceConfiguration()
configuration(::CurrentVolume) = CurrentConfiguration()
configuration(::ReferenceBoundary) = ReferenceConfiguration()
configuration(::CurrentBoundary) = CurrentConfiguration()

same_configuration(a::AbstractConfiguration, b::AbstractConfiguration) = typeof(a) === typeof(b)
same_configuration(a::AbstractMeasure, b::AbstractConfiguration) = same_configuration(configuration(a), b)
same_configuration(a::AbstractConfiguration, b::AbstractMeasure) = same_configuration(a, configuration(b))
same_configuration(a::AbstractMeasure, b::AbstractMeasure) = same_configuration(configuration(a), configuration(b))
