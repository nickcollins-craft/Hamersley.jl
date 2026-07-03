abstract type AbstractTrace end

struct ValueTrace <: AbstractTrace end
struct NormalGradientTrace <: AbstractTrace end
struct NormalComponentTrace <: AbstractTrace end
struct TangentialComponentTrace <: AbstractTrace end
struct EdgeValueTrace <: AbstractTrace end

requires_higher_continuity(::AbstractTrace) = false
requires_higher_continuity(::NormalGradientTrace) = true
