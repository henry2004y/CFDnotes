@def title = "Mesh Generation"
@def hascode = true
@def date = Date(2020, 6, 16)
@def rss = "Mesh generation in AE623 notes."
@def tags = ["syntax", "code"]

\toc

# Unstructured Mesh Generation

Interface: [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl)

[Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl)

##

[TriangleMesh.jl](https://github.com/konsim83/TriangleMesh.jl): Delaunay and constraint Delaunay meshes.

[VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl): 2D Delaunay and Voronoi tessellations on generic point types.

# Structured Mesh Generation

Generally, a mapping is required from the reference domain to the physical domain.

## Multi-block

I haven't found any implementation yet in Julia.
