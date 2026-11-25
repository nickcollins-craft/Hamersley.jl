#=
The Hamersley test suite.

Code-quality checks (Aqua and JET) come first, then the structural tests, from the core
vocabulary out to linearity and coupling. The shared test problems are defined once, in
`test_problems.jl`.
=#


using Hamersley
using Ferrite: Ferrite
using Test
using Aqua
using JET


include("test_problems.jl")


@testset "Hamersley.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Hamersley)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(Hamersley; target_modules = (Hamersley,))
    end
    include("core_tests.jl")
    include("expression_tests.jl")
    include("material_interface_tests.jl")
    include("formulation_tests.jl")
    include("trace_tests.jl")
    include("interpolation_tests.jl")
    include("linearity_coupling_tests.jl")
    include("integrator_tests.jl")
end
