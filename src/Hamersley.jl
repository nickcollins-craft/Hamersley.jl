module Hamersley

include("core/configurations.jl")
include("core/traces.jl")
include("core/fields.jl")
include("core/operators.jl")
include("core/quantities.jl")

include("responses/responses.jl")
include("responses/stress.jl")
include("responses/flux.jl")
include("responses/storage.jl")
include("responses/source.jl")

include("forms/weak_terms.jl")
include("forms/boundary_conditions.jl")
include("forms/equations.jl")
include("forms/problems.jl")
include("forms/audit.jl")

include("backends/ferrite.jl")
include("prelude.jl")

end
