#=
Generating test geometries.

Generators for typical geomechanics test geometries (cylindrical triaxial specimens,
oedometer cells, strip footings, …). Each one builds the geometry with the Gmsh API,
meshes it, names the physical groups that boundary conditions will refer to, and returns a
`Ferrite.Grid`. The contract is fixed here; the generators themselves are not yet
written.
=#


"""
    cylindrical_specimen(; radius, height, element_size, order = 1)

Generate a cylindrical specimen, as used in triaxial tests, of the given `radius` and
`height` (m), meshed with elements of characteristic size `element_size` (m) and geometric
order `order`.

# Returns
A `Ferrite.Grid` with facet sets `:top`, `:bottom` and `:lateral`.
"""
function cylindrical_specimen(;
    radius::Real,
    height::Real,
    element_size::Real,
    order::Integer = 1,
)
    not_yet_implemented("generation of cylindrical specimens")
    return nothing
end
