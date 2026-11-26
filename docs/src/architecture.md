# Architecture

This page records how Hamersley is organised, and tries to give some explanation as to _why_
it is organised that way.

Hamersley is a finite element framework for higher-order and multiphysics solid mechanics,
built on [Ferrite](https://ferrite-fem.github.io/Ferrite.jl/stable/). It is designed to be
able to handle finite strains, and solves them in a total Lagrangian manner.

## Layers

The source is organised in layers, each of which only uses the ones above it.

| Directory         | Responsibility                                                                     |
| ----------------- | ---------------------------------------------------------------------------------- |
| `core/`           | Configurations, kinematic assumptions, quantity tags, traces                       |
| `materials/`      | The constitutive law interface and the material library                            |
| `expressions/`    | The dataflow graph: fields, operators, kinematics, transforms, constituive laws    |
| `formulation/`    | Weak terms, equations, boundary conditions, problems, the formulation verification |
| `discretisation/` | Interpolations, grids, dofs, quadrature                                            |
| `state/`, `time/` | Quadrature-point state; time integrators and rate factors                          |
| `compile/`        | The single step from a formulation to type-stable element kernels                  |
| `assembly/`       | Assembly modes, element kernels, global assembly                                   |
| `solvers/`        | Newton and cached-factorisation strategies                                         |
| `io/`             | Output                                                                             |
| `meshing/`        | The `Meshing` submodule, the only place that loads Gmsh                            |

The initial layers (`core` to `formulation`) describe the continuous problem and check that
the formulation makes sense. The compile step turns a checked formulation into explicitly
typed kernels so that the solver can proceed as fast as possible.

## Kinematic assumptions

Every problem needs to declare `kinematics = SmallStrain()` or
`kinematics = FiniteStrain(motion = :u)`. Here, `motion = :u` indicates that the
displacement field is named `:u`, and it is being used for the finite strain. If the
displacement field has a different name (e.g. `:d`), that should be used instead. Equally,
it can be worthwhile to name fields to clarify their role (e.g. in a Biot-type formulation,
one may use `:u_s` to make clear that it is the solid motion that defines the
configuration).  
Any operation that depends on the distinction between configurations but isn't explicitly
specified will result in the verification throwing an error.

- Under `SmallStrain`, `F ≈ I` and `J ≈ 1`: the reference and current configurations
  coincide, all of the stress measures also coincide with the Cauchy stress, and the
  pull-back is simply the identity (so it doesn't do anything). Cauchy stress paired with
  `ε` over the reference volume is valid and will pass the verification check.
- Under `FiniteStrain`, the reference and current configurations are distinct from each
  other, and every integral is taken over the reference configuration, requiring pull-backs
  to bring the terms into that configuration.

## Formulations as a dataflow graph

A formulation (when looked at from a very high level) is a small directed graph:

```
fields → operators → tagged quantities → constitutive laws → tagged quantities → weak terms
```

The fields are the quantities that are interpolated by the finite elements (displacement,
temperature etc), operators are applied to those fields (gradients, divergence etc), tagged
quantities are the stresses and strains (and their analogues for other fields) which are
linked by constitutive laws (there is a tagged quantity that is input into the constitutive
law, and a tagged quantity that is an output of the consitutive law), and the weak terms are
the componenets of the overall formulation that the finite element method discretises.  
The **tag** in a tagged quantity says what it _represents_ (which stress or strain measure,
which configuration it is measured per unit of, how it transforms). There are many different
possibilities for where its value comes from. For example, it may be interpolated directly
from a field (`Value(:σ)` for a field declared with `represents = CauchyStress()`), derived
by an operator (`DeformationGradient(:u)`), or returned by a constitutive law
(`output(binding, :P)`). The formulation verification checks conjugacy and configurations
based on the tags alone, so mixed formulations don't require any special treatment.

Each weak term names what elements are playing the role of test functions explicitly with
`Variation`: `Term(Variation(DeformationGradient(:u)), P, ReferenceVolume())` is
`∫ P : δF \, \mathrm{d}v_0`, and
`Term(Variation(GreenLagrangeStrain(:u)), S, ReferenceVolume())` is
`∫ S : δE \, \mathrm{d}v_0`. They are different terms with different linearisations.

## Linearisation

- **Linearisation needs to pass through the whole graph, not only the constitutive laws.**
  Each vertex supplies its derivative with respect to the fields, and the compile step
  applies the chain rule from each weak term back to the fields. Typically these links are
  trivial: field operators are linear, and under small strains so is every kinematic vertex,
  leaving the constitutive laws as the only nonlinearity (in the case that they are
  nonlinear, they could be linear too!). Under finite strains the linearisation of the
  kinematic vertices is not trivial: the geometric stiffness comes from `δE` (equivalently
  `F` in `P = F S`) depending on `u`, and the load stiffness from pulled-back follower
  loads.
- **Laws typically differentiate their time-discrete update.** For history-dependent
  materials the tangent is (in general) the algorithmic tangent, with respect to the
  current-step inputs. A local return mapping provides its Jacobian so that the implicit
  function theorem gives the consistent tangent. However, it is possible (and occasionally
  desirable) to write a different tangent within the constitutive law. There are
  circumstances where a theoretically sub-optimal tangent being returned by the constitutive
  law will result in desirable global behaviours. This is left to the expertise of the user.
- **Laws declare dependencies per output.** The global sparsity pattern and the off-diagonal
  blocks of a multiphysics tangent are built from these, so structurally zero blocks can be
  declared. This helps build efficient assembly routines.
- **Rates take their factor from the time integrator.** A `Rate` vertex is an ordinary
  constitutive law input; the integrator supplies the relevant time integration factor.
- **What can be differentiated is separate from what a solver actually uses.** A staggered
  scheme may use only the diagonal blocks whereas a classical Newton--Raphson solve will use
  them all.
- **Tangent symmetry is declared.** A constitutive law says whether the tangent it returns
  is symmetric (`tangent_is_symmetric`). This is useful as it potentially allows the
  exploitation of symmetry in the global solver. However, a symmetric constitutive tangent
  does not by itself make the global system symmetric: follower loads and non-symmetric
  multiphysics couplings can still break the symmetry of the global stiffness matrix and
  thus preclude exploiting the symmetric system solvers.

## Pull-backs and push-forwards

There is a `PullBack` vertex and a `PushForward` vertex. The specific transformation rule
that is applied depends on the quantity's tag, so users don't need to write it. Under small
strains both push-forward and pull-back vertices are the identity.

The formulation verification rejects a quantity integrated over the wrong configuration and
suggests the transform, but it won't insert one without any prompting: where a transform
goes (before or after a constitutive law, on a quantity or on its rate) is a modelling
decision. Under finite strains, rates are taken of reference quantities (`\dot{E}`,
`\dot{J}`, `\dot{T}`), which are objective automatically.

## Boundary traces

The traces of each field are derived from the operators acting on its test function in the
volume terms: a gradient leaves the value, a second gradient leaves the value, the normal
gradient and edge values, and a divergence leaves the normal component. Natural conditions
infer their trace from the quantity they prescribe. Whether an essential condition can be
imposed strongly follows from Ferrite's conformity of the interpolation. The current
implementation only supports Lagrangian elements, which support a strong imposition, but
this is not the case for every element available in Ferrite.

Unspecified traces default to homogeneous natural conditions (that can be thought of as
analogous to "if a boundary doesn't specify a traction or a displacement, it experiences a
zero traction"). `explain(problem)` reports the status of every trace on every named
boundary, and `strict_boundaries = true` turns the default non-standard trace into a
formulation verification error, so that explicitly accounting for every boundary can be
imposed if need be. Formulation templates (for example a future second-gradient template
discretised micromorphically) translate continuum boundary conditions onto the discretised
fields.

## Elements

Fields declare a Ferrite family and an order (`Field(:u, VectorField(), Lagrange, 2)`). The
family is resolved against each subdomain's reference shape at compile time. Currently the
system admits only Lagrange elements, but extending to admitting another family should only
need a new method of `is_supported_family`, because:

- conformity is read from Ferrite rather than declared;
- only scalar families are vectorised;
- kernels evaluate through Ferrite's cell values, reinitialise with the cell itself, and
  never treat degrees of freedom as nodal values (in assembly, boundary conditions, output,
  or mass lumping).

Element-local fields (discontinuous pressures or stresses) are part of the structure for
static condensation, which is not yet implemented.

## Linear problems

Linearity is inferred (`is_affine`) or claimed (`linear = true`); a claim is checked once,
by the formulation verification, and then trusted. A linear problem uses the same backend as
a nonlinear problem, as it is equivalent to one (guaranteed-to-converge) Newton step with a
cached factorisation, which is exact because the residual is affine. The factorisation is
reused for as long as its key (time step, parameters, constrained set) is unchanged;
changing essential values only changes the right-hand side of the problem, so factorisations
can be reused without issue in this case. Kernels can assemble the residual only, the
tangent only, or both.

## Dependencies

- Ferrite, Tensors, LinearAlgebra and SparseArrays are core dependencies.
- Gmsh and FerriteGmsh are ordinary dependencies, used only by the `Meshing` submodule,
  whose sole output is a `Ferrite.Grid` with named sets. Meshes from anywhere else are just
  as good.
- NonlinearSolve is not a dependency. Material-point solves go through `solve_local`; the
  choice of solver behind it will be made when plasticity models are written.
