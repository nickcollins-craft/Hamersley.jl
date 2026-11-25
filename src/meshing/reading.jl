#=
Reading meshes.

We read Gmsh files through FerriteGmsh, which turns Gmsh's physical groups into Ferrite's
named cell and facet sets.
=#


"""
    read_mesh(path)

Read the Gmsh mesh at `path` and return it as a `Ferrite.Grid`, with each physical group
of the mesh available as a named cell or facet set.
"""
function read_mesh(path::AbstractString)
    # We check for the file ourselves, because Gmsh's own error is not informative.
    if !isfile(path)
        throw(ArgumentError("no mesh file found at $path"))
    end
    grid = FerriteGmsh.togrid(String(path))
    return grid
end
