using Hamersley
using Test
using Aqua
using JET

@testset "Hamersley.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(Hamersley)
    end
    @testset "Code linting (JET.jl)" begin
        JET.test_package(Hamersley; target_defined_modules = true)
    end
    # Write your tests here.
end
