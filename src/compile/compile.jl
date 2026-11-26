#=
Compilation.

The formulation layer is deliberately not type-stable: it is small, it is built once, and
it is easier to read and to check as a loose graph. The hot loop of assembly must be
type-stable. Compilation is the single step between the two: it verifies the problem,
checks it against the grid, resolves every interpolation per subdomain, and turns the
graph into concretely typed element kernels (closures over Ferrite cell values, with
fields addressed through named tuples).

This replaces the earlier "Ferrite backend" abstraction. Hamersley is built on Ferrite,
so Ferrite's types are used directly below this point.
=#


"""
    CompiledProblem

A problem compiled against a grid: degree-of-freedom handler, constraint handler,
concrete element kernels, sparsity pattern and storage for internal state. The fields of
this type will be fixed when compilation is written.
"""
struct CompiledProblem end


"""
    compile(problem, grid)

Compile `problem` against `grid`.

We first verify the formulation of the problem (throwing on any issue) and check that the grid provides every
set the problem names; both checks are implemented. Building the kernels is not yet.
"""
function compile(problem::Problem, grid::Ferrite.AbstractGrid)
    # The formulation must be sound before anything is discretised.
    verification(problem; throw_on_error=true)

    # The grid must provide every set the problem names.
    validate_grid(grid, problem)

    not_yet_implemented("compilation of a problem into element kernels")
    return CompiledProblem()
end
