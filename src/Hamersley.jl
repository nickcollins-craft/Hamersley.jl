#=
Hamersley: finite elements for higher-order and multiphysics solid mechanics, built on
Ferrite.

The package is organised in layers, each of which only uses the ones above it:

- `core/`: the vocabulary (configurations, kinematic assumptions, quantity tags, traces,
  regions and prescribed values);
- `materials/`: the constitutive law interface and the material library;
- `expressions/`: the dataflow graph of a formulation (fields, operators, kinematics,
  pull-backs, constitutive bindings);
- `formulation/`: weak terms, equations, boundary conditions, problems and the formulation verification;
- `discretisation/`: interpolations, grids, degrees of freedom and quadrature;
- `state/` and `time/`: internal state and time integration;
- `compile/`: the single step from a formulation to type-stable element kernels;
- `assembly/` and `solvers/`: element kernels, global assembly and solution strategies;
- `io/`: output;
- `meshing/`: the `Meshing` submodule, the only place that loads Gmsh.

The architecture is described in full in `docs/src/architecture.md`.
=#


module Hamersley

using Ferrite: Ferrite, Lagrange
using LinearAlgebra
using SparseArrays
using Tensors

# The vocabulary: errors, configurations, kinematic assumptions, quantity tags, traces,
# regions and prescribed values.
include("core/errors.jl")
include("core/configurations.jl")
include("core/kinematic_assumptions.jl")
include("core/field_shapes.jl")
include("core/quantities.jl")
include("core/traces.jl")
include("core/regions.jl")
include("core/prescriptions.jl")

# Constitutive laws only refer to quantity tags, so they come before the graph.
include("materials/interface.jl")
include("materials/placeholder.jl")
include("materials/local_solve.jl")
include("materials/library/isotropic_linear_elasticity.jl")
include("materials/library/neo_hookean.jl")

# The dataflow graph of a formulation.
include("expressions/fields.jl")
include("expressions/vertices.jl")
include("expressions/kinematics.jl")
include("expressions/transforms.jl")
include("expressions/constitutive.jl")
include("expressions/boundary_vertices.jl")
include("expressions/traversal.jl")

# Weak forms, problems and the formulation verification.
include("formulation/terms.jl")
include("formulation/equations.jl")
include("formulation/boundary_conditions.jl")
include("formulation/schedules.jl")
include("formulation/templates.jl")
include("formulation/problems.jl")
include("formulation/condition_terms.jl")
include("formulation/trace_derivation.jl")
include("formulation/verification.jl")

# Discretisation.
include("discretisation/interpolations.jl")
include("discretisation/grids.jl")
include("discretisation/dofs.jl")
include("discretisation/quadrature.jl")

# Internal state and time integration.
include("state/quadrature_state.jl")
include("time/integrators.jl")

# Compilation, assembly and solution.
include("compile/linearity.jl")
include("compile/coupling.jl")
include("compile/compile.jl")
include("assembly/kernels.jl")
include("assembly/assemble.jl")
include("solvers/newton.jl")
include("solvers/linear.jl")
include("io/output.jl")

# Meshing is a submodule, so that Gmsh is only ever loaded inside it.
include("meshing/Meshing.jl")
using .Meshing: read_mesh

include("exports.jl")

end
