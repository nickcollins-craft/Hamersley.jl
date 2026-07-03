using Hamersley
using Test
using Aqua
using JET

@testset "Hamersley.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Hamersley)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(Hamersley; target_modules=(Hamersley,))
    end
    include("core_tests.jl")
    include("response_tests.jl")
    include("boundary_condition_tests.jl")
    include("audit_tests.jl")
end
