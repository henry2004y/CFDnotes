@def title = "Mesh Generation"
@def hascode = false
@def date = Date(2020, 6, 16)
@def rss = "Mesh generation in AE623 notes."
@def tags = ["syntax", "code"]

\toc

# Unstructured Mesh Generation

Interface: [MeshIO.jl](https://github.com/JuliaIO/MeshIO.jl)

[Gmsh.jl](https://github.com/JuliaFEM/Gmsh.jl)

About 2D Delaunay triangulation:
* brute force algorithm has $\mathcal{O}(n^2)$ complexity, but the explanation in the note is wrong. See a better explanation [here](https://cs.stackexchange.com/questions/2400/brute-force-delaunay-triangulation-algorithm-complexity).

##

[TriangleMesh.jl](https://github.com/konsim83/TriangleMesh.jl): Delaunay and constraint Delaunay meshes.

[VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl): 2D Delaunay and Voronoi tessellations on generic point types.

# Structured Mesh Generation

Generally, a mapping is required from the reference domain to the physical domain.

## Multi-block

I haven't found any implementation yet in Julia.

I noticed the so called H grid, which can produce high-quality elements near the object surface. Why are we not using it?

## Mappings

Algebraic Map: Figure 2.4.13 presents a highly concise C code snippet for transfinite interpolation. So elegant.

PDE-based mapping: still a mystery to me.

# Libraries

By far the best library regarding mesh in Julia:
[Meshes.jl](https://github.com/JuliaGeometry/Meshes.jl)

Built-in mesh support (including tree-based AMR):
[Trixi.jl](https://github.com/trixi-framework/Trixi.jl)