#=
Internal state.

Internal variables (plastic strains, damage, hardening variables) live at quadrature
points, not in finite element fields. They are stored twice: the converged state of the
previous step (`old`) and the state of the current iteration (`new`), so that a failed
step can be rolled back and a converged one committed.

Element-local unknowns that have been statically condensed out of the global system are
stored the same way, per element rather than per quadrature point.
=#


"""
    QuadratureState(old, new)

Internal state at every quadrature point of every cell.

# Fields
- `old::Vector{Vector{S}}`: the converged state of the previous step, indexed as
  `old[cell_index][quad_point_index]`.
- `new::Vector{Vector{S}}`: the state of the current iteration, indexed likewise.
"""
struct QuadratureState{S}
    old::Vector{Vector{S}}
    new::Vector{Vector{S}}
end


"""
    commit!(state)

Overwrite `state.old` with `state.new`, once a step has converged.
"""
function commit!(state::QuadratureState)
    not_yet_implemented("committing quadrature-point state")
    return nothing
end


"""
    rollback!(state)

Overwrite `state.new` with `state.old`, when a step is cut back.
"""
function rollback!(state::QuadratureState)
    not_yet_implemented("rolling back quadrature-point state")
    return nothing
end


"""
    ElementLocalState(old, new)

Statically condensed element-local unknowns, recovered element by element after each
global solve.

# Fields
- `old::Vector{Vector{T}}`: the converged values of the previous step, per element.
- `new::Vector{Vector{T}}`: the values of the current iteration, per element.
"""
struct ElementLocalState{T}
    old::Vector{Vector{T}}
    new::Vector{Vector{T}}
end
