#=
The meshing submodule.

Mesh reading and the generation of typical geomechanics test geometries live here, and
only here: this is the one place in Hamersley that loads Gmsh. Everything this submodule
produces is a `Ferrite.Grid` with named cell and facet sets, which is the only thing the
rest of the package accepts. A mesh from anywhere else (Ferrite's `generate_grid`, a
user's own Gmsh script, another mesher) is therefore just as good, and the submodule could
later move to a sibling package without changing anything else.
=#


"""
    Hamersley.Meshing

Mesh reading and generation, returning Ferrite grids with named cell and facet sets. This
is the only part of Hamersley that depends on Gmsh.
"""
module Meshing

using Ferrite: Ferrite
using FerriteGmsh: FerriteGmsh
using Gmsh: Gmsh

import ..Hamersley: not_yet_implemented

include("reading.jl")
include("generation.jl")

export read_mesh, cylindrical_specimen

end
