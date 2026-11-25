# Architecture

This page records how Hamersley is organised, and why. It is the design record: when a
decision here changes, this page changes with it.

Hamersley is a finite element framework for higher-order and multiphysics solid mechanics,
built on [Ferrite](https://ferrite-fem.github.io/Ferrite.jl/stable/). It starts from small
strains but is designed from the outset for finite strains, solved in a total Lagrangian
manner.

## Layers

The source is organised in layers, each of which only uses the ones above it.

| Directory         | Responsibility                                                       |
|-------------------|----------------------------------------------------------------------|
| `core/`           | Configurations, kinematic assumptions, quantity tags, traces         |
| `materials/`      | The constitutive law interface and the material library              |
| `expressions/`    | The dataflow graph: fields, operators, kinematics, transforms, laws  |
| `formulation/`    | Weak terms, equations, boundary conditions, problems, the audit      |
| `discretisation/` | Interpolations (and the element gate), grids, dofs, quadrature       |
| `state/`, `time/` | Quadrature-point state; time integrators and rate factors            |
| `compile/`        | The single step from a formulation to type-stable element kernels    |
| `assembly/`       | Assembly modes, element kernels, global assembly                     |
| `solvers/`        | Newton and cached-factorisation strategies                           |
| `io/`             | Output                                                               |
| `meshing/`        | The `Meshing` submodule, the only place that loads Gmsh              |

The formulation layers (`core` to `formulation`) describe and check a problem; they are
not type-stable and do not need to be. The compile step turns a checked formulation into
concretely typed kernels for the hot loop, so type stability is paid for once, in one
place.

## Kinematic assumptions

Every problem declares `kinematics = SmallStrain()` or `kinematics = FiniteStrain(motion =
:u)`. Nothing that depends on the distinction between configurations is ever guessed.

- Under `SmallStrain`, ``F ≈ I`` and ``J ≈ 1``: the reference and current configurations
  coincide, every stress measure collapses onto the Cauchy stress, and every pull-back is
  the identity. Cauchy stress paired with ``ε`` over the reference volume is valid.
- Under `FiniteStrain`, configurations are distinct, and every integral is taken over the
  reference configuration. "Total Lagrangian" is how finite strains are solved, not the
  name of the assumption.

## Formulations as a dataflow graph

A formulation is a small directed graph:

```
fields → operators → tagged quantities → laws → tagged quantities → weak terms
```

Each quantity carries a **tag** saying what it *means* (which stress or strain measure,
which configuration it is measured per unit of, how it transforms). Where its value comes
from is a separate matter: it may be interpolated from a field (`Value(:σ)` for a field
declared with `represents = CauchyStress()`), derived by an operator
(`DeformationGradient(:u)`), or returned by a law (`output(binding, :P)`). The audit
checks conjugacy and configurations on tags alone, so mixed formulations need no special
treatment.

Each weak term names its test side explicitly with `Variation`:
`Term(Variation(DeformationGradient(:u)), P, ReferenceVolume())` is ``∫ P : δF \,
\mathrm{d}V_0``, and `Term(Variation(GreenLagrangeStrain(:u)), S, ReferenceVolume())` is
``∫ S : δE \, \mathrm{d}V_0``. They are different terms with different linearisations.

## Linearisation

- **Every node has a derivative**, not only laws. Under finite strains the geometric
  stiffness comes from kinematic nodes (``P = F S``, ``δE`` depending on ``u``), and load
  stiffness from pulled-back follower loads.
- **Laws differentiate their time-discrete update.** For history-dependent materials the
  tangent is the algorithmic tangent, with respect to the current-step inputs only. A local
  return mapping provides its Jacobian so that the implicit function theorem gives the
  consistent tangent.
- **Laws declare dependencies per output.** The global sparsity pattern and the
  off-diagonal blocks of a multiphysics tangent are built from these, so structurally zero
  blocks stay out of the pattern.
- **Rates take their factor from the time integrator.** A `Rate` node is an ordinary law
  input; the integrator supplies ``a = ∂\dot{x}/∂x_{n+1}`` (``1/Δt`` for backward Euler,
  ``γ/(βΔt)`` and ``1/(βΔt^2)`` for Newmark). Laws never contain it.
- **What can be differentiated is separate from what a solver uses.** A staggered scheme
  uses only the diagonal blocks; full Newton uses them all.
- **Potentials are declared.** A law deriving from ``Ψ`` has a symmetric tangent.

## Pull-backs and push-forwards

There is one `PullBack` node and one `PushForward` node. The transformation rule (Piola for
fluxes, ``P = Jσ F^{-T}`` for stresses, covariant for gradients, ``J`` for volume
densities, Nanson's formula for area densities) is a property of the quantity's tag, so
users never name it. Under small strains both nodes are the identity, so a formulation
written with them runs under either assumption.

The audit rejects a quantity integrated over the wrong configuration and suggests the
transform, but never inserts one: where a transform goes (before or after a law, on a
quantity or on its rate) is a modelling decision. Under finite strains, rates are taken of
reference quantities (``\dot{E}``, ``\dot{J}``, ``\dot{θ}``), which are objective
automatically.

## Boundary traces

The traces of each field are derived from the operators acting on its test function in
the volume terms: a gradient leaves the value, a second gradient leaves the value, the
normal gradient and edge values, and a divergence leaves the normal component. Natural
conditions infer their trace from the quantity they prescribe. Whether an essential
condition can be imposed strongly follows from Ferrite's conformity of the interpolation;
weak (Nitsche) imposition is not yet implemented.

Unspecified traces default to homogeneous natural conditions. `explain(problem)` reports
the status of every trace on every named boundary, and `strict_boundaries = true` turns a
defaulted non-standard trace into an audit issue. Formulation templates (for example a
future second-gradient template discretised micromorphically) translate continuum
boundary conditions onto the discretised fields.

## Elements

Fields declare a Ferrite family and an order (`Field(:u, VectorField(), Lagrange, 2)`).
The family is resolved against each subdomain's reference shape at compile time, behind a
single gate that currently admits only Lagrange. Admitting another family should only
need a new method of `is_supported_family`, because:

- conformity is read from Ferrite rather than declared;
- only scalar families are vectorised;
- kernels evaluate through Ferrite's cell values, reinitialise with the cell itself, and
  never treat degrees of freedom as nodal values (in assembly, boundary conditions,
  output, or mass lumping).

Element-local fields (discontinuous pressures or stresses) are part of the structure for
static condensation, which is not yet implemented.

## Linear problems

Linearity is inferred (`is_affine`) or claimed (`linear = true`); a claim is checked once,
by the audit, and then trusted. A linear problem is not a separate assembly path: it is
one Newton step with a cached factorisation, exact because the residual is affine. The
factorisation is reused for as long as its key (time step, parameters, constrained set) is
unchanged; changing Dirichlet values only changes the right-hand side. Kernels can
assemble the residual only, the tangent only, or both.

## Dependencies

- Ferrite, Tensors, LinearAlgebra and SparseArrays are core dependencies.
- Gmsh and FerriteGmsh are ordinary dependencies, used only by the `Meshing` submodule,
  whose sole output is a `Ferrite.Grid` with named sets. Meshes from anywhere else are just
  as good.
- NonlinearSolve is not a dependency. Material-point solves go through `solve_local`; the
  choice of solver behind it will be made when plasticity models are written.
