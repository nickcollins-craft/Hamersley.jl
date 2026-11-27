# Boundary conditions: design

!!! note "Status"
    This page is a design. Steps 1 to 5 of the [Order of implementation](@ref) (the
    formulation layer) are implemented; steps 6 to 11 need the runtime and are planned.
    The [Status](@ref) section says exactly what exists, and the sections below mark what
    is still planned.

This page sets out how Hamersley should specify boundary conditions at the level of
sophistication of the earlier `finite_element_system` code (magnitudes and directions that
vary in time, component and point selection, periodicity, and conditions that change type
part-way through an analysis), without giving up what the formulation layer already does
well: deriving which traces need conditions, checking them before anything is
discretised, and knowing which loads follow the deformation. It also covers what the
earlier code could not do: switching when a measured quantity reaches a threshold, follower
pressures, and rigid platens that rotate with the specimen.

## Status

**Implemented** (the formulation layer, with tests):

- **Regions**: `FacetSet`, `EdgeSet` and `VertexSet`; a bare symbol is a facet set.
- **Component selections**: `components = 2` or `components = (1, 3)` on every condition,
  with a status per component and inheritance per component.
- **Value prescriptions**: `Constant`, `Prescribed(f; breakpoints, current_position)` and
  `Tabulated(times, values)`. `PlaceholderValue` is gone, and `NaturalBC` takes a `value`.
- **Follower loads**: the `ReferencePressure` and `SpatialPressure` tags; natural conditions
  compiled to boundary terms (over the current boundary for a follower load), so that
  `is_follower_load` is retired; values of follower loads that depend on the current
  position.
- **Robin conditions**: `RobinBC`, compiled to a boundary term; boundary terms written in
  the equations count in trace status.
- **The conflict rule**, per phase and per component.
- **Schedules**: `Phase`, `Schedule`, `AtTime`, `Never` and `Release`, with inheritance by
  key, verification per phase, `is_affine(problem, phase)` and
  `default_strategy(problem, phase)`, and the coupling pattern taken over every phase.
- `explain` reports the trace statuses (and the discrete conditions) phase by phase.

**Planned** (needs the runtime): `Hold()`, `PrescribedRate`, global unknowns and rigid
bodies, `PeriodicBC`, `When` triggers and the step controller, `Normal()` and
`Tangential()`, line loads and concentrated forces, checking the shapes of `Prescribed`
values at compile time, and marking the tangent non-symmetric for follower loads.

**Where the implementation differs from the first design**, the sections below have been
updated. In short:

- Ferrite grids carry no named edge sets, so an edge region is usually given as the
  meeting of two facet sets, `EdgeSet(:top, :front)`.
- Trace statuses are reported as `Prescribed`, `Loaded`, `LoadedByTerm` and
  `DefaultedToZero` (the types behind them are `PrescribedTrace`, `LoadedTrace`,
  `TermLoadedTrace` and `DefaultedToZero`, since `Prescribed` is now a value prescription).
- Strict mode reports a non-standard trace only when no condition covers any of its
  components; whether a remainder such as "all components except 3" is empty depends on
  the spatial dimension, so that check belongs to the compile step.

## What the earlier code does

`finite_element_system` attaches to each condition a value callback `(x, ctx)`, which may
also be rate- or acceleration-controlled (integrated with the time scheme's stage
weights); a component selection (Dirichlet only); a region (`:facet` or `:vertex`);
periodic pairs of facet sets with an optional transform (homogeneous only); and an
activation window `active(t)`. A condition whose window opens can latch the degrees of
freedom of the committed solution, which is how "load by traction, then hold the
displacement reached" is written. The switch itself happens at a prescribed time:
switching when a load or displacement reaches a threshold is not supported.

The value callbacks receive the committed solution and material states through `ctx`.
That is powerful, but it cannot be checked ahead of time, and a value that depends on the
solution makes the residual depend on the solution in a way the tangent never sees.

## Principles

1. **Separate *what*, *where*, *how much* and *when*.** The condition types (*what*) stay as
   they are. A region and a component selection (*where*), a value prescription (*how
   much*) and a phase schedule (*when*) are separate objects that conditions are built
   from.
2. **Values only ever affect the right-hand side**, with one principled exception. A
   prescribed value enters the residual (natural conditions) or the constrained values
   (essential conditions), and never the tangent, so no value prescription can change
   whether a problem is linear. The exception is a follower load, whose dependence on the
   deformation is declared by its tag and linearised with it (see
   [Follower loads](@ref)). Loads that depend on the solution in any other way (convection,
   Robin conditions, elastic supports) are boundary *terms* of the weak form, where the
   linearity, coupling and measure analyses see them.
3. **Declare the schedule; don't hide it in callbacks.** Opaque `active(t)` functions
   cannot be verified. A list of phases with explicit triggers can: every phase is
   checked as a problem in its own right before the analysis starts.
4. **Reading the solution is explicit and rare.** The one way a prescribed value reads the
   solution is `Hold()`, evaluated once when its phase starts and constant afterwards.
5. **Constraints are exact, and assumptions are checked rather than enforced.** Hamersley
   uses neither penalties nor a contact solver. Constraints are imposed by eliminating
   degrees of freedom (strongly, or through affine constraints). Where a model rests on an
   assumption that a contact solver would otherwise enforce (a platen stays in contact and
   does not slip), the assumption is checked after each step and a violation is reported.
6. **The runtime state of the boundary is stored like material state**: an old and a new
   copy, committed when a step converges and rolled back when it is cut.

## *Where*: regions and components

`on` accepts a region. A bare symbol keeps meaning a facet set, so existing problems are
unchanged.

| Region                   | Meaning                                         | Ferrite set |
| ------------------------ | ----------------------------------------------- | ----------- |
| `FacetSet(:top)`, `:top` | Faces in 3D, edges in 2D (codimension one)      | facet set   |
| `EdgeSet(:top, :front)`  | Edges in 3D where two facet sets meet           | two facet sets |
| `EdgeSet(:crease)`       | A named edge set in 3D                          | none yet    |
| `VertexSet(:corner)`     | Individual points, for example to pin rigid modes | vertex set |

`EdgeSet` only exists in 3D, where it is also where the edge forces of second-gradient
continua act (`EdgeForce`, conjugate to `EdgeValueTrace`); in 2D the facets are already
edges, and the grid check reports an `EdgeSet` as an error. Ferrite grids carry facet and
vertex sets but no named edge sets, so the grid check reports a named `EdgeSet` too, and
asks for the meeting of two facet sets instead; named edge sets wait for the mesh readers
to provide them.

Facet conditions are the ones trace derivation is about. Edge and point conditions on an
ordinary continuum are discrete rather than variational: they never appear when the weak
form is integrated by parts, so they do not change the status of any facet trace. Their
essential conditions are imposed strongly (Ferrite's `Dirichlet` accepts `EdgeIndex` and
`VertexIndex` sets). Their natural conditions need their own tags: a line load (force per
unit length) on an edge, and a concentrated force at a point, both conjugate to
`ValueTrace`. These are planned; until then the verification reports an area load on a
vertex or edge region.

Components are selected with `components`, which applies to the *trace*, not to the field:
for `ValueTrace` of a vector field they are components of ``u``, and for
`NormalGradientTrace` components of ``∂_{n_0} u``.

```julia
EssentialBC(:u; on = :bottom, components = (3,), value = 0.0)
```

Two further selections, planned, are for sliding on vector ``H^1`` fields: `Normal()` (prescribe
``u \cdot N``) and `Tangential()` (prescribe the tangential part). Both use the
**reference** normal ``N``. On a flat face this is exact for any amount of deformation and
sliding, because a rigid plane's normal does not change: ``N \cdot (X + u) = c`` is
linear in ``u``. On faces aligned with the axes these are component constraints; on
oblique faces they are one `AffineConstraint` per node, reported as a third imposition,
`AffineImposition()`. Constraints that follow the *current* orientation of a face are the
subject of [Rigid bodies and rigid ties](@ref).

**Status per component.** Trace status is a status per (region, field, trace,
component). This is what lets one face be displacement-controlled in ``z`` and
traction-loaded in ``x`` and ``y``, and it is what the conflict rule below checks.
Because the number of components is only known with the mesh, selections are stored
symbolically (every component, some components, or every component except some), and
what no condition covers is reported the same way:

```
right: u ValueTrace() [components 1] — Prescribed
right: u ValueTrace() [all components except 1] — DefaultedToZero
```

## *How much*: value prescriptions

A value prescription replaces `PlaceholderValue`, and `NaturalBC` gains a `value`
keyword. Its meaning follows from the tag: a `ReferenceTraction` value is a force per unit
reference area, a `SpatialTraction` value a force per unit current area.

| Prescription                                       | Meaning                                       |
| -------------------------------------------------- | --------------------------------------------- |
| `Constant(v)`                                      | ``v``, everywhere on the region and at all times |
| `Prescribed(f; breakpoints = (), current_position = false)` | ``f(X, t)``, or ``f(X, x, t)`` (see [Follower loads](@ref)) |
| `Tabulated(times, values)`                         | piecewise linear in time, uniform in space    |
| `PrescribedRate(f; order = 1, initial)` (planned)  | ``\dot{v} = f(X, t)``                         |
| `PrescribedRate(f; order = 2, initial, initial_rate)` (planned) | ``\ddot{v} = f(X, t)``        |
| `Hold()` (planned)                                 | the value at the start of the phase (see below) |

Plain numbers and functions are accepted as shorthand for `Constant` and `Prescribed`.
``X`` is the reference position, as everywhere in a total Lagrangian formulation.

**Why `Tabulated` and `breakpoints`.** Any of these could be written as `Prescribed(f)`;
the others exist because the system can see inside them. The step controller lands on
every breakpoint, as it does on `AtTime` triggers, because a kink in the loading history
smeared across a step loses accuracy for rate-dependent materials. A `Tabulated` value
declares its breakpoints (its rows) and its shape (uniform in space) automatically; a ramp
is a two-row table, and cyclic or square-wave histories are longer ones. A `Prescribed`
function declares its kinks through `breakpoints`.

**Prescribed rates.** `PrescribedRate` prescribes a time derivative of order `order` of
whatever the condition prescribes: a displacement, a temperature, a normal gradient or a
load. The name and the `order` keyword follow the `Rate(vertex; order)` vertex and the
`RateOf` tag. The prescription is integrated by the problem's own time integrator, through
the same update that relates a `Rate` vertex to the unknowns, so that a rate on the
boundary and a rate in the bulk are integrated consistently. The formulation verification
checks that the integrator defines the order (`BackwardEuler` does not define order 2).
`initial_rate` is required for order 2 and an error for order 1, and both `initial` and
`initial_rate` accept `Hold()`.

**Shapes.** The shape of a value must match the trace and the selected components. For
`Constant` and `Tabulated` this is checked by the formulation verification (a scalar for a
scalar trace or a pressure; a tuple for a vector trace or a traction, of the length of the
selection when components are selected; a single selected component may be given as a
number). A function is to be checked at compile time, by evaluating it once on the region.
Tables are also checked for rows that match their times, and times that strictly
increase. `prescribed_value(prescription, X, t)` evaluates any prescription, and
`breakpoints(prescription)` returns the times a step controller should land on.

### `Hold()` (planned)

`Hold()` makes a key continuous across a change of phase. It reads the committed solution
once, when its phase starts, from a converged state, and only for the key it sits on; for
the rest of the phase it is a constant, so it adds nothing to the tangent, cannot change
linearity or coupling, and rolls back trivially. What it reads follows from the
work-conjugate pair on the key, and from what held the key in the previous phase:

| Previous phase on this key | `Hold()` on an essential condition | `Hold()` on a natural condition |
| -------------------------- | ---------------------------------- | ------------------------------- |
| essential                  | the trace value it was held at     | the reaction (the quantity conjugate to the trace) |
| natural                    | the trace value reached            | the boundary quantity it was carrying |

The trace value may be a displacement, a temperature, a pore pressure or a normal gradient;
the conjugate quantity a force, a heat flow, a fluid flux or a double force.

For a dead load, holding the reaction as the consistent nodal vector and holding it as a
field are the same thing. For a follower load they are not: a frozen nodal vector fixes the
force at each node in magnitude and direction, which turns the follower load into a dead
load. So a held follower load is held as a field. When the previous phase was natural the
field is the prescription itself, evaluated at the start of the phase. When it was
essential, the reaction is projected onto the boundary with the consistent boundary mass in
the current configuration (``∫ N_a N_b \, \mathrm{d}a \; t_b = r_a``), which reproduces
the reaction exactly at the switch and then follows the deformation. Holding a pressure
keeps only the normal part of the reaction; the runtime reports the tangential part it
discards when it is not negligible.

## Follower loads

A load is a follower when it is given per unit current area, or when its direction follows
the deformation. Three changes make these first-class.

1. **Tags.** `SpatialTraction` stays for a force per unit current area in a fixed spatial
   direction. A pressure, the most common follower load, gets its own tag, `SpatialPressure`:
   its value is a scalar ``p``, its traction ``-p \, n`` per unit current area, and it is
   conjugate to `ValueTrace`. (`BoundaryPressure` is already taken by mixed ``H(\mathrm{div})``
   formulations, where it is conjugate to `NormalComponentTrace`.)
2. **Natural conditions are compiled to boundary terms.** A `NaturalBC` becomes
   `Term(Variation(Value(field)), load, measure)`: over `ReferenceBoundary` for a dead
   load, and over `CurrentBoundary` (Nanson's formula, ``n \, \mathrm{d}a = J F^{-T} N \,
   \mathrm{d}A``) for a follower load. The load stiffness then comes from the same
   linearisation as every other term, and the measure-aware linearity and coupling analyses
   treat it without a special case, so `is_follower_load` is retired.
3. **Values may depend on the current position.** A follower load's value may depend on the
   current position ``x = X + u`` as well as on ``(X, t)``, as a hydrostatic pressure
   ``p = ρ g (h - x_3)`` on a deforming body does. The dependence is declared and
   differentiable, so it enters the load stiffness. Dead loads, and every essential value,
   remain functions of ``(X, t)`` only.

A pressure's load stiffness is not symmetric in general (it is for closed surfaces and
some support conditions), so a follower load is to mark the global tangent as
non-symmetric, so that a symmetric factorisation is not chosen; this is planned with the
assembly.

Natural conditions on facet regions are compiled (`condition_term`); an edge force of a
second-gradient continuum is not, because there is no edge measure yet, but it is a dead
load and changes neither linearity nor coupling. The test side of a term whose trace is not
the value (a double force, a pressure on a flux field) is a `BoundaryTrace(field, trace)`. Under `SmallStrain` the follower tags reduce
to their reference counterparts through `canonical_quantity`, and the problem is linear
again.

## Robin conditions and other boundary terms

Conditions that depend on the solution throughout a phase are boundary terms of the weak
form. The common one is a boundary condition of its own,

```julia
RobinBC(:T; on = :outer, coefficient = h, reference = T_∞, measure = ReferenceBoundary())
```

compiled into the term ``∫ δT \, h (T - T_∞) \, \mathrm{d}A`` (over the current boundary
with `measure = CurrentBoundary()`, when ``h`` is given per unit current area). For a vector
field it is an elastic support, ``∫ δu \cdot k (u - u_∞) \, \mathrm{d}A``. The coefficient
is a constant, since it multiplies the unknown and so enters the tangent; the reference
value is any prescription, since it only enters the right-hand side. The term is affine. Trace status counts any boundary term whose test side is
`Variation(Value(field))` as loading that field's value trace on its facet set, so the
face is no longer reported as defaulted to zero, and it is not an error in strict mode. A
boundary term may share a key with a natural condition (a convective face with an imposed
heat flux, say), since both are boundary integrals and their sum is explicit in the weak
form; a boundary term on a key that an essential condition prescribes is reported, and so
are two Robin conditions on one key. A boundary term over the whole boundary (no `on`)
loads whatever no condition covers, and is never in conflict, since writing ``∫_{∂Ω}``
while part of the boundary is constrained is ordinary practice.

## *When*: phases and triggers

A schedule is a sequence of phases. Each phase lists the conditions that change and the
trigger that ends it.

```julia
Phase(name; conditions = (...), until = AtTime(t) | When(observable, threshold) | Never())
```

**Inheritance and replacement.** Each phase inherits the conditions of the previous one. A
condition replaces an inherited one with the same key (region, field, trace, component),
and `Release(field; on, components)` returns a key to the homogeneous natural default.
Swapping a face from traction to displacement control is then a single condition in the
new phase, and because keys are per component, an essential condition on ``u_z``
replaces only the ``z`` component of an inherited traction. `explain` prints every phase
in full, so that the inherited set never has to be reconstructed by hand.

**Triggers.**

- `AtTime(t)`: the step controller (planned) shortens the last step of the phase so that it lands
  exactly on ``t``.
- `When(observable, threshold)` (planned): the phase ends when a monitored quantity crosses a value.
  Observables are functionals declared in the formulation: the resultant reaction on a
  region (`Reaction(:u; on = :bottom, component = 3)`), the average of a trace over a
  region (`Average(:u; on = :top, component = 3)`), the resultant of a natural condition,
  or a reaction or position of a rigid body. When a step crosses the threshold, the step
  controller locates the crossing by cutting the step (bisection or a secant on ``Δt``) to
  within a tolerance.
- `Never()`: the last phase runs to the end of the analysis.

**What happens at a phase boundary.** At the start of a phase the runtime resolves every
`Hold()` from the committed solution, rebuilds the constraint handler if the constrained
set has changed, and continues. Because `FactorisationKey` already includes the set of
constrained degrees of freedom, a cached factorisation is refreshed automatically.

**One phase is the default.** `Problem(...; boundary_conditions = (...))` stays as it is,
and means a single phase, called `default`, that never ends.
`Problem(...; schedule = Schedule(phases...))` is the general form; giving both is an
error. `phase_conditions(problem, phase)` returns the conditions of a phase (by name or
position) after inheritance.

**Linearity and strategy per phase.** `is_affine(problem, phase)` and
`default_strategy(problem, phase)` look at one phase; `is_affine(problem)` and
`default_strategy(problem)` require every phase to be affine. The coupling pattern is the
union over every phase, so that it is built once.

## Rigid bodies and rigid ties

This section is planned; it needs [Global unknowns](@ref).

A constraint that follows the *current* orientation of a face, such as "no tangential
displacement in the plane of the rotated face", cannot be written node by node.
Displacement is measured from the reference position, so a rigid rotation of the face
already gives every node a tangential displacement; and the current tangent plane belongs
to the deformed surface, so it depends on the neighbouring nodes. What such a condition
means physically is that the face does not slip relative to something rigid: a platen.

A platen that the face stays in contact with and sticks to moves the face rigidly. So
under those two assumptions, contact against a rigid platen is exactly a rigid tie, and
needs neither a contact solver nor a penalty.

**Rigid bodies.** `RigidBody(:top_platen)` declares a body with its own unknowns: a
translation ``u_c`` and a rotation ``θ`` (three and three in 3D, two and one in 2D). These
are global unknowns, not cell fields (see [Global unknowns](@ref)).

**Rigid ties.** `RigidTie(:u; on = :top, body = :top_platen)` ties every point of the
region to the body,

```math
u(X) = u_c + \left( R(θ) - I \right) (X - X_c),
```

with ``X_c`` the centroid of the region unless given. The tie is imposed by elimination:
the face degrees of freedom are expressed in terms of the body's. Under small rotations
``R ≈ I + [θ]_×`` and the tie is affine (a set of `AffineConstraint`s with fixed
coefficients), so under `SmallStrain` the problem stays linear. Under finite rotations the
tie is nonlinear; its linearisation,
``δu = δu_c + δθ × R(θ)(X - X_c)``, refreshes the constraint coefficients at every Newton
iterate, and the second derivative of ``R``, contracted with the face's reactions, adds a
rotational geometric stiffness for quadratic convergence. For trace status, the tie
prescribes the value trace of the region; it couples the face's degrees of freedom to the
body's in the sparsity pattern.

**Conditions on a body.** Each unknown of a body takes an ordinary essential or natural
condition, with the same value prescriptions, phases and `Hold()` as everything else, and
is free (natural with zero force or moment) unless one is given. A `frame` keyword says
whether components are in the laboratory frame (`Lab()`, a piston pushing along a fixed
axis) or in the body's own rotated frame (`Body()`, where the translation is
``u_c = R(θ) \tilde{u}`` and conditions apply to ``\tilde{u}``). A force along the body's
own axis is a follower force on the body, with its load stiffness through ``θ``.

| Laboratory condition                 | Condition on the body                               |
| ------------------------------------ | --------------------------------------------------- |
| displacement control along the axis  | essential on the axial translation                  |
| force control along the axis         | natural on the axial translation (one resultant force) |
| ball seat (free to tilt)             | natural, zero moment, on the tilting rotations (the default) |
| fixed seat                           | essential on the rotations                          |
| laterally guided platen              | essential on the tangential translations            |

Because these are single unknowns, `Hold()` on a body holds one reaction or position, the
follower complications of holding a field do not arise, and `When` triggers can name the
quantities measured in the laboratory directly.

**Assumption checks.** Contact and stick are assumptions, so they are checked after every
converged step rather than enforced: `RigidTie(...; check = StickAndContact(μ = 0.5))`
recovers the traction on the face from the reactions, and reports where and when the
normal traction stops being compressive (lift-off) or the tangential traction exceeds
``μ`` times the normal traction (slip). The model stays free of contact, and it says when
it stops being valid.

**Limits.** This models rigid platens. A flexible membrane, slip once friction is
exceeded, or contact between deformable bodies are contact problems, which Hamersley does
not solve. A face that slides along a curved rigid surface (``φ(X + u) = 0``) is a
nonlinear equality constraint that needs multipliers on the boundary, and so waits for
boundary fields (see [Known structural gaps](@ref)).

### Example: a triaxial test with a rotating top platen

```julia
platen = RigidBody(:top_platen)

schedule = Schedule(
    Phase(
        :isotropic_compression;
        conditions = (
            EssentialBC(:u; on = :bottom, value = (0.0, 0.0, 0.0)),
            RigidTie(:u; on = :top, body = :top_platen, check = StickAndContact(μ = 0.5)),
            # The cell pressure acts on the deformed membrane, so it is a follower load.
            NaturalBC(
                :u;
                on = :lateral,
                quantity = SpatialPressure(),
                value = Tabulated((0.0, 1.0), (0.0, 200.0e3)),
            ),
            NaturalBC(
                :top_platen;
                component = Translation(3),
                frame = Body(),
                value = Tabulated((0.0, 1.0), (0.0, -200.0e3*specimen_area)),
            ),
        ),
        until = AtTime(1.0),
    ),
    Phase(
        :shearing;
        # Switch the platen from force to displacement control where it stopped. The
        # rotations stay free (a ball seat), and the cell pressure is inherited.
        conditions = (
            EssentialBC(
                :top_platen;
                component = Translation(3),
                frame = Body(),
                value = PrescribedRate((X, t) -> -1.0e-4; initial = Hold()),
            ),
        ),
        until = When(Position(:top_platen, Translation(3); frame = Body()), -0.02),
    ),
)

problem = Problem(;
    fields = (u,),
    equations = (momentum,),
    rigid_bodies = (platen,),
    schedule = schedule,
    kinematics = FiniteStrain(motion = :u),
)
```

## Periodicity

This section is planned.

Periodicity constrains two regions against each other, so it is neither essential nor
natural on either one.

```julia
PeriodicBC(field; mirror, image, components = all, transform = nothing, offset = nothing)
```

- `transform` maps mirror points to image points; by default it is found from the facet
  sets (Ferrite's `collect_periodic_facets`).
- `offset` makes the condition inhomogeneous:
  ``u(X_\text{image}) = u(X_\text{mirror}) + \bar{G}\,(X_\text{image} - X_\text{mirror})``,
  where ``\bar{G}`` is a macroscopic gradient. Under `SmallStrain` it is ``\bar{ε}``;
  under `FiniteStrain`, total Lagrangian, it is ``\bar{F} - I`` in reference coordinates.
  A prescribed gradient is given as a value prescription (`MacroGradient(value)`), and
  Ferrite's `PeriodicDirichlet` accepts it for periodicity aligned with the axes; other
  cases go through `AffineConstraint`s. A gradient that is partly unknown
  (stress-controlled homogenisation) makes its unknown components global unknowns, with
  the average stress as their conjugate conditions.
- Both regions get a new trace status, `Paired()`, for the paired components: their
  boundary terms cancel, so neither is reported as defaulted to zero, and in strict mode
  neither is an error.
- Pinning a point of a periodic boundary (the usual way to remove rigid translations) is
  allowed; the compile step applies it to the master degree of freedom of the pair.
- The constraints couple degrees of freedom across the body, so the sparsity pattern
  includes them (Ferrite's `add_constraint_entries!`).

## Global unknowns

This section is planned.

Rigid bodies and unknown macroscopic gradients both need unknowns that belong to no cell:
a handful of scalars, each with an ordinary essential or natural condition, coupled to the
cell degrees of freedom through constraints. Ferrite's `DofHandler` and
`ConstraintHandler` only know about cell degrees of freedom, so Hamersley owns a thin
layer that appends the global unknowns to the system, builds the elimination of the tied
or periodic degrees of freedom in terms of them, and adds their entries to the sparsity
pattern. Both features are built on it, so it is written once.

## Verification

The formulation verification checks every phase, after inheritance, as a problem in its
own right; when there are several phases, a message that belongs to one names it ("in
phase shearing, …"). Every condition is checked once, as it was written. Rule 8 now covers
the following.

**Implemented:**

1. **Conflicts.** In each phase, every key (region, field, trace, component) has at most
   one claim that excludes the others: an essential condition together with anything, two
   natural conditions, or two Robin conditions on overlapping components are reported.
   (Two tractions on one key could be superposed, but an explicit sum is clearer than a
   silent one.) Discrete conditions (two pins of the same vertex set) are checked too.
2. **Schedules.** Phase names are unique; `AtTime` triggers strictly increase; only the
   last phase may never end; a `Release` must remove something the phase inherits.
3. **Regions.** A condition acts on a kind of region that can carry it: area loads and
   Robin conditions on facets, edge forces on edges, a value anywhere, other traces on
   facets. Values at points and along edges are discrete, and are not checked against the
   derived traces (so a pressure field can be pinned at a point).
4. **Components and values.** Selections are well formed, scalars are not given
   components, and constant and tabulated values have the right shape.
5. **Follower values.** Only a follower load may depend on the current position.
6. **Robin conditions.** The field has a Robin flux (a scalar or vector field), and the
   measure is a boundary measure.

**Planned:** that the time integrator defines the order of every `PrescribedRate`; that a
`Hold()` in the first phase has an initial solution to read; the rules on rigid bodies
(every tied body is declared, every condition on a body names one of its unknowns,
`frame = Body()` only for translations and forces).

## Known structural gaps

Two pieces of structure are missing, and each blocks more than one feature.

- **Global unknowns**, needed for rigid bodies and unknown macroscopic gradients. This page
  plans them.
- **Fields that live only on the boundary**, needed for multipliers on the boundary: a face
  sliding on a curved rigid surface, and any future contact or mortar coupling. This page
  does not plan them.

## Effects on the rest of Hamersley

| Part                        | Change                                                     | Status   |
| --------------------------- | ---------------------------------------------------------- | -------- |
| trace status, rule 8        | per phase, per component; counts boundary terms and Robin conditions; edge and vertex regions | done |
| `is_follower_load`          | retired: natural conditions are compiled to boundary terms, and the measure-aware analyses see them | done |
| `grid_issues`               | checks vertex sets and edge regions as well as facet sets  | done     |
| `is_affine`                 | a problem is affine only if every phase is; per phase on request | done |
| `default_strategy`          | chosen per phase on request                                | done     |
| `coupling_matrix`           | the union over all phases                                  | done     |
| trace status                | counts rigid ties and periodic pairs                       | planned  |
| `default_strategy`          | follower loads make the tangent non-symmetric              | planned  |
| sparsity                    | constraint entries and the global unknowns                 | planned  |
| `state/`                    | a `BoundaryState` with `old` and `new` copies (rate-integrated values, held values, the current phase) | planned |
| `time/`                     | a step controller that lands on `AtTime` triggers and breakpoints, and locates `When` crossings | planned |
| Meshing                     | named edge sets, and helpers that add the standard bounding-box facet, edge and corner sets to a grid | planned |

None of this changes how a weak form is written. The formulation layer has gained new
types, and the static analyses loop over phases; the runtime is to gain global unknowns, a
boundary state and a step controller.

## Order of implementation

Formulation-layer work, now done:

1. The conflict rule, and boundary terms (Robin conditions) in trace status, for the single
   phase that existed before schedules; `RobinBC`.
2. Regions (`FacetSet`, `EdgeSet`, `VertexSet`) and component selection, with status per
   component.
3. Value prescriptions (`Constant`, `Prescribed`, `Tabulated`), replacing
   `PlaceholderValue`, and `value` on `NaturalBC`.
4. Natural conditions compiled to boundary terms, `SpatialPressure`, values that depend on
   the current position for follower loads, and `is_follower_load` retired.
5. Schedules with `AtTime` triggers, inheritance and `Release`, verified per phase; the
   union over phases for coupling and linearity.

Work that needs the runtime (assembly, the constraint handler, the step controller), still
to do:

6. `Hold()`, with the boundary state.
7. `PrescribedRate`, through the time integrators.
8. Global unknowns, then rigid bodies, rigid ties, frames and assumption checks.
9. `PeriodicBC`, with prescribed and then unknown offsets.
10. `When` triggers, with observables and event location.
11. `Normal()` and `Tangential()` on oblique faces; line loads and concentrated forces.

## Decisions

These were agreed in discussion and are recorded here so that they are not reopened by
accident.

- Phases **inherit** conditions, and replace them by key.
- Prescribed values see ``(X, t)`` only, except for follower loads (which may also see the
  current position) and `Hold()`.
- Two natural conditions on one key are **not superposed**; a boundary term may share a key
  with a natural condition.
- The solver strategy is chosen **per phase**.
- **No penalties and no contact solver.** Constraints are exact, by elimination; the
  assumptions that replace contact are checked after each step.
