# Hamersley.jl style guide

## The guiding principle

**The code should read as much like a paper as possible.** Someone with the paper (or
thesis, or derivation notes) open beside the source should be able to follow each line. To
aid this, there are three broad principles:

1. **Use Mathematical notation.** Names, symbols and the "visual shape" of expressions match
   the written maths.
2. **Lots of prose.** Docstrings and comments explain what is being computed and why, in
   full sentences.
3. **Performance-based choices should be explained.** High-performance representations
   (preallocation, in-place updates, packed arrays) are used deliberately where worthwhile,
   and are documented so they don't cost readability.

When a rule otherwise conflicts with that principle, use the principle. If not specified in
this guide, try and follow the surrounding code.

---

## 1. Language

1.1. **Write in Australian English** in all prose: docstrings, comments, error messages,
documentation and commit messages.

1.2. Identifiers follow the same spelling where they are English words (`initialise_system`,
`normalised_residual`). The exception is names fixed by an external API, which keep the
upstream spelling (`LinearAlgebra.normalize` is called as it is).

1.3. Write prose in the first-person plural, as in a paper: "We first predict the elastic
step", "Here we protect against…".

---

## 2. Naming

2.1. **Use the paper's symbols, in Unicode.** Variables that correspond to a quantity in the
maths take its symbol: `γ_e`, `κ_e`, `Ψ`, `σ_c`, `T`, `Δε_v`, `η_L`.

2.2. **Subscripts become underscore-separated parts**: `p_t`, `δ_c_n`, `γ_t`.

2.3. **Various additional accent marks, indices etc can be written as suffixes**, e.g.:

| Suffix                   | Meaning                                         | Example     |
| ------------------------ | ----------------------------------------------- | ----------- |
| `_dot`                   | time rate                                       | `γ_t_dot`   |
| `_bar`                   | barred quantity                                 | `M_hat_bar` |
| `_hat`                   | hatted quantity                                 | `i_hat`     |
| `_T`                     | transpose                                       | `H_n_T`     |
| `_inv`                   | inverse (or a factorisation standing in for it) | `M_hat_inv` |
| `_nought`                | initial / reference value                       | `φ_nought`  |
| `_acc`                   | accumulated over time                           | `γ_p_acc`   |
| `_old`, `_new`, `_trial` | time-step levels                                | `u_old`     |

2.4. **Mathematical capitalisation beats language convention.** Matrices and tensors are
capitals when they are capitals in the paper (`K`, `M_hat`, `H_n`, `S_inv`), and
mathematical objects named after people keep their proper names for clarity (`Kronecker_δ`,
`Heaviside`, `Macaulay`).

2.5. **Everything that is not a symbol is long, descriptive snake_case English**:
`current_state_variables`, `extended_stress_control_vector`, `facet_quad_point_index`.
Abbreviate only with well-established and unambiguous domain terms (`czm`, `LCP`, `dof`,
`quad`, `bc`).

2.6. **Loop variables are descriptive** (`time_step`, `element_index`, `czm_counter`). The
exception is index notation, where `i, j, k, l` are chosen to mirror the mathematical
convention.

2.7. **Role suffixes** can be used to make a variable's job obvious:

- `_val`: the value of the function of the same name (`Ψ_val = Ψ(state)`)
- `_vals`, `_history`: storage over time steps
- `_helper`: an intermediate term in a long piece of algebra (`A_helper`)
- `_mask`: a Boolean selector
- `_local`, `_cache`: working copies and preallocated buffers

<!-- End of role suffixes list. -->

2.8. **Types and modules are `CamelCase`** (`FirstPiolaStress`, `EssentialBC`), following
Julia convention. Spell them out in English rather than with symbols.

2.9. **Public API names are ASCII.** Exported functions, types and keyword arguments are
things users must type, so give them ASCII names (`yield_function`, `time_step = h`) to
facilitate that. Unicode symbols should be reserved for internals usage such as local
variables, struct fields, and non-exported helpers.

2.10. **Mutating functions end in `!`**, and the mutated arguments come first:
`stress!(kinetic_variables, dΨ_val)` mutates kinetic_variables (it may or may not mutate
dΨ_val).

---

## 3. Docstrings

Every function, type and constant that is more than a trivial accessor gets a docstring.
Docstrings say _what_ the object is and _where it comes from_; comments inside the body say
_how_ and _why_ (see section 4).

3.1. **Layout.** Use the standard Julia form: the signature indented by four spaces, a
one-sentence summary, then whatever else is needed.

````julia
"""
    yield_function(kinetic_variables, state_variables)

Evaluate the yield surface ``y`` for the current stress and internal state.

The yield surface couples the state locator ``F`` with the deviatoric stress invariant
``q_e``, normalised by the critical-state stress ratio ``M``:

```math
y = F^2 + \\left( \\frac{q_e}{M p_e + φ (1 - B) c} \\right)^2 - 1
```

A negative value means the trial state is elastic.

# Arguments
- `kinetic_variables`: packed kinetic variables (layout in [`KineticLayout`](@ref)).
- `state_variables`: packed state variables ``[γ_e, κ_e, ρ, φ, B, λ]``.

# Returns
The scalar value of ``y``.

# References
<Authors> (<year>), Eq. (<n>).
"""
function yield_function(kinetic_variables, state_variables)
    # …
end
````

3.2. **Write the maths in the docstring.** Use ` ```math ` blocks and double-backtick inline
maths so Documenter renders them. If the docstring has a lot of LaTeX, use
`@doc raw""" … """` rather than escaping every backslash.

3.3. **Cite the source where possible.** Say which paper, chapter or equation the function
implements (`# References`), so the reader can put the derivation beside the code.

3.4. **Mutating functions say what they mutate**, and which entries: "Overwrites entries
22:30 of `kinetic_variables` with the total stress."

3.5. **Struct docstrings describe every field**: its symbol, its meaning, and its units or
dimension where relevant.

3.6. **Files open with a short block comment** (`#= … =#`) describing what the file contains
and how it fits into the package.

---

## 4. Comments

4.1. **Comment densely.** As a rough guide, about one line in four or five should be a
comment, and almost every one to three lines of code should sit under one. A body reads as a
series of short paragraphs, each led by a comment that says what that paragraph does.

4.2. **Comments are full sentences**, in Australian English and the first-person plural.

4.3. **Give the reason, not only the action.** Put the reason in parentheses or a second
clause: "(we use a small numerical value for consistency with the libplast code)", "−1.0
because the normal points outwards but `H_n` needs the inward normal".

4.4. **Write the equation next to the code that computes it** when the code is not obviously
the equation: `# ∂y/∂γ·γ_p_bar`.

4.5. **Label defensive code** so it isn't mistaken for physical reasoning e.g. "Protect
against nonsensical values of B arising from round-off".

4.6. **Document layouts where they are used.** When a packed array is filled or read, say
which entries hold what (or use named index ranges, see 8.4, which makes this
self-documenting).

4.7. **Section banners** split long files into coherent sections to facilitate organisation
e.g.:

```julia
# ---------------------------------------------------------------------------------------- #
# ----------------------- Functions that implement the return map ------------------------ #
# ---------------------------------------------------------------------------------------- #
```

---

## 5. Layout and whitespace

5.1. **Lines are at most 92 characters**, including comments and docstrings. The only
exceptions are things that cannot be broken, such as URLs.

5.2. **Break long signatures one argument per line**, with a trailing comma:

```julia
function implicit_return_map(
    state::PlasticityState,
    caches::ReturnMapCaches,
    loading::LoadingProtocol,
)
    # …
end
```

5.3. **Break long maths at term boundaries (if needed), inside parentheses, with the
operator leading each continuation line**, so it reads like a displayed equation:

```julia
residual_φ = (
    φ_new*(1 - (1 - χ)*Δε_v)
    - φ_current
    - h*λ*φ_p_bar
)
```

The parentheses are required. Outside them, a line that ends with a complete expression ends
the statement, and `- φ_current` on the next line becomes a separate, silently discarded
statement.

5.4. **Operator spacing follows written maths.** Multiplication, division and powers don't
have spacing around them. Addition, subtraction, assignment and comparisons are spaced:

```julia
i_hat = M*v_old - h*θ*K*u_old + h*θ*(θ*F + (1 - θ)*F_old)
```

so `*` is visually much tighter than `+`, as it does on paper. This facilitates reading the
terms as an equation.

5.5. **Keyword arguments and defaults have spaces around `=`**:
`solve(prob; reltol = 1e-12)`, `function Problem(; equations::Tuple = ())`.

5.6. **Blank lines**: two between top-level definitions, one between paragraphs inside a
function, and one either side of a section banner.

---

## 6. Expressions

6.1. **Use brackets generously** where they help to make grouping obvious: `((h*θ)^2)`,
`(1/3)*(σ_11 + σ_22 + σ_33)`.

6.2. **Avoid augmented assignment for scalars**: write `t = t + h` and `k = k + 1`, not
`t += h`. This again helps the interpretation of the material.

6.3. **For arrays, keep the explicit form but broadcast in place**: `u .= u .+ h .* v` (or
`@. u = u + h*v`), not `u = u + h*v`, which allocates a new array.

6.4. **Negation and literals are explicit in concrete code**: `-1.0*normal`, comparisons
against `0.0`.

6.5. **In generic code, literals respect the element type.** Use `zero(T)`, `one(T)`,
`zeros(T, n)` and `similar(x)` rather than `0.0` and `zeros(Float64, n)`, so the code works
with `Float32` and automatic-differentiation dual numbers.

6.6. **Name the element type when allocating concrete storage**: `zeros(Float64, (9, n))`,
`falses(18)`.

---

## 7. Control flow

7.1. **Spell out piecewise definitions as `if`/`elseif`/`else`**:

```julia
# Macaulay bracket
function Macaulay(x)
    if x < zero(x)
        return zero(x)
    else
        return x
    end
end
```

A library function (`clamp`, `max`) is fine if it is the clearer statement, but don't squash
more complicated or physically meaningful case splits into a one-liner.

7.2. **Returns are explicit.** Every function that returns a value ends in `return`, often
through a named result (`Ψ_val = …; return Ψ_val`). In-place functions end in
`return nothing`.

7.3. **Choose behaviour by type, not by string.** In the package, prefer dispatch on types
(`ReferenceConfiguration`, `DrainedTriaxial()`) over string flags
(`load_type == "drained_triaxial"`). Use `Symbol`s for names and labels.

---

## 8. Data and state

8.1. **Group related data in structs.** A function that needs many related values takes a
struct holding them, not a long positional argument list. As a rule of thumb, more than five
or six positional arguments is a sign that some belong together.

8.2. **Make structs concretely typed** (parametric where needed) so they stay type-stable:

```julia
struct MaterialParameters{T<:Real}
    K_bar::T  # Bulk modulus (Pa)
    G_bar::T  # Shear modulus (Pa)
    ζ::T      # Coupling between breakage and porosity
end
```

8.3. **Choose the state representation to fit the problem.** Flat preallocated vectors with
a fixed layout are the right choice for tight loops: constitutive updates at a Gauss point,
nonlinear solves, automatic differentiation through a residual. They are not the default for
everything. Use structs, `Tensors.jl`/`StaticArrays` types or named tuples wherever they are
clearer and don't cost measurable performance.

8.4. **When you do pack state into a flat vector, name the index ranges**, so the code
doesn't depend on magic numbers that only comments explain:

```julia
const ELASTIC_STRESS = 1:9
const COUPLE_STRESS = 10:18
const YIELD_SURFACE = 60

kinetic_variables[ELASTIC_STRESS] = dΨ_val[ELASTIC_STRESS]
```

Document the full layout once, in a docstring or a layout struct.

8.5. **No non-`const` globals.** Parameters arrive through arguments (usually a parameter
struct). Global values that really are fixed are `const`.

8.6. **Copying is explicit.** Use `copy`/`deepcopy` when keeping an old value, and `x .= y`
(or `copyto!`) for in-place overwrites. Never rely on aliasing by accident.

---

## 9. Performance

9.1. **Make high-performance choices deliberately, and comment them.** Preallocate in hot
loops, update in place, and use AD-compatible caches (`PreallocationTools.DiffCache`,
`similar(x, …)`). Say in a comment why a choice was made when it costs readability.

9.2. **Lean on established libraries for numerics**: ForwardDiff, NonlinearSolve, NLopt,
Ferrite, SparseArrays and so on. Don't rewrite a solver that a maintained package already
provides.

9.3. **State tolerances explicitly** at the call site (`reltol = 1e-12, abstol = 1e-12`).

---

## 10. Errors and numerical defences

10.1. **Always throw.** Use `throw(ArgumentError("…"))`, `throw(DomainError(…))` or
`error("…")`. Never construct an exception without throwing it: `ErrorException("…")` on its
own does nothing.

10.2. **Error messages list the valid options** when an input is not recognised: "Invalid
loading type, please choose from: biaxial, drained_triaxial, …".

10.3. **Check solver return codes** (`SciMLBase.successful_retcode`, NLopt return symbols)
and fail loudly, reporting the code.

10.4. **Test for NaN with `isnan`**, never `x == NaN`, which is always `false`.

10.5. **Defend against round-off where the physics demands it**: clamp quantities that must
stay in a physical range, and fall back to finite differences when an AD gradient comes back
non-finite. Label every such defence (4.5).

---

## 11. Duplication

11.1. **A block of code that appears twice should be a function.** Small, single-purpose
functions named after the quantity they compute are the norm. Copy-pasted blocks drift apart
when one copy is fixed and the other is not.

11.2. **An orchestrating function calls the small ones in a documented order**, with a
comment for each step, so it reads like an algorithm described by a paper.
