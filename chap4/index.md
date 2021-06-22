@def title = "The Discontinuous Galerkin Method"
@def hascode = true
@def date = Date(2020, 6, 16)
@def rss = "Discontinuous Galerkin method in AE623 notes."
@def tags = ["syntax", "code"]

\toc

# One-Dimensional Conservation Laws

## 1D Scalar Conservation

The full Julia implementation can be found at ...
Note that when p ≥ 4, the scheme becomes unstable. I wonder what's the reason behind.

## 1D Advection-Diffusion

# Libraries

[MultivariatePolynomials.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl): package for constructing Legendre basis functions. It is not very user-friendly.

[FEMQuad.jl](https://github.com/JuliaFEM/FEMQuad.jl): numerical integration schemes for cartesian and tetrahedron domains. I use it for obtaining Gauss-Legendre integral points in triangles, i.e. Dunavant points.

# Code Optimizations

Even if some programming languages have libraries for polynomials and basis, it is generally not good to evaluate polynomials everytime when you do computations within the time steps. Remember the extremely slow version of 1D shallow-water DG solver implemented with symbolic functions in MATLAB? I tend to write that way the first time because it has a very similar form to the mathematical expression, but in practice, I should not do that.

For example, let the coefficients be stored in an array named `coef_PVE`. To get the density at a point location with basis function `ϕ` in reference space `(ξ, η)`,
```
ρ = sum(coef_PVE[:,ρ_,  k] .* ϕ)(ξ=>x[1], η=>x[2])
```

Alternatively, we can precompute the basis function values at quadrature points and store them as float arrays. Since the unknowns of FEM are multiplier coefficients of basis functions, we can have highly efficient computations of getting physical quantities and evaluating integrals through simple floating point arithmatics:
```
ρ = sum(coef_PVE[:,ρ_,   k] .* phi[:,j])
```
where `phi` is the precomputed values of basis functions evaluated at quadrature points.

However, since the solutions of finite element methods are functions, in principle we can obtain values at any given location in our simulation domain. When doing visualizations, it may be better to evaluate our solution at more points than the quadrature points used for computation. In this case, we can switch back to the polynomial evaluations. 