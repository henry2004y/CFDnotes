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

Originally I had implemented my own RK4 timestepping method. Inspired by Trixi.jl and the idea of semi-discretization, I managed to use DifferentialEquations.jl directly for timestepping, which is much more robust since it is not limited to one time marching scheme. This is a much more elegant solution, especially for numerical packages.

## 1D Advection-Diffusion


# Two-Dimensional Conservation Laws

Key points:
* The Lagrangian basis points in the reference triangle is not the same as quadrature points.
* The area integral quadrature points in the reference triangle are not the same as line integral quadrature points along the triangle edges. See the [discussion](https://scicomp.stackexchange.com/questions/27441/line-integral-along-the-edge-of-an-isoparametrically-mapped-triangle). The integration direction also matters: by convention, counter-clockwise path integral is positive while clockwise path integral is negative. One easy check you can do is to integrate a uniform function along the edge and compare that with the perimeter of the triangle.
* Line integral mapping is error-prone compared with area integral, because you have to take direction and edge lengths into consideration. In area integral from reference element to physical element, you just need to multiply by the determinant of the Jacobian matrix; in line integral, you need to multiply by the line length as well as be careful about whether the quadrature points should go counter-clockwise or clockwise. Generally, we need counter-clockwise quadrature points for the left element of the edge evaluation, and clockwise quadrature points for the right element of the edge evaluation. Maybe there are tricks taking advantage of symmetry since the reference domain in 1D is usually [-1,1], but I have not had a working version with that yet.
* Pay attention to the necessary order of quadrature rules. In my original Fortran code, when evaluating the edge fluxes I use Gauss-Legendre integral of order `2p+1`, which may be unnecessarily high given that the evaluated polynomial basis functions do not reach that order. If you happen to use an inapppropriate low order quadrature rule (i.e. Dunuvant points and weights), you may end up with a singular elemental mass matrix!

For the first time, it took me a month to get working version with Fortran. I generally followed the structure design of BATSRUS, which was a good learning process. I thought at the time that it was clean, well-written and well-documented, until four years later when I tried to rewrite the model in Julia. For the second time, it took me about a week to get a working version with Julia. I spent quite some time reading my own Fortran code (not so easy, even I coded the entire thing myself~), and found many places that were either wrong (e.g. L2 error computation) or could be improved (e.g. directory checking/creation, variable naming, module structures). Then it took me a day to optimize the code with 10-20x speed up by removing small temporary array allocations.
The profiling results are summarized [here](https://github.com/henry2004y/CFD/issues/1).

I do not finish everything in the project:
* least square projection for initialization is there, but never being tested. However, it is more important to understand why this works and how it works. This leads to the question of how to evaluate the accuracy of finite element solutions. Instead of being more accurate on a certain bunch of points, we may want to be more accurate on the whole domain instead. The least square projection follows the idea that if I generate the solution by minimizing the L2 norm error, then it should provide more accurate solution later on when evaluating the L2 norm error.
* Plotting is quite inaccurate in the Julia implementation because for simplicity I only have one averaged value per element. Unlike in Matlab we have this `fill/patch` method which can draw per element so we can clearly observe the discontinuities, in Matplotlib we only have `tricontour` or `tripcolor`. `tricontour` only accepts nodal values; `tripcolor` can accept either cell values or node values. The simplest way to get a relatively accurate plot is probably averge the solution on each node and use gradient fill for the elements. It would be even better if I can find a way to plot element-by-element.

## Libraries

[MultivariatePolynomials.jl](https://github.com/JuliaAlgebra/MultivariatePolynomials.jl): package for constructing Legendre basis functions. It is not very user-friendly.

[TypedPolynomials.jl](https://github.com/JuliaAlgebra/TypedPolynomials.jl): for evaluating numerical values of multivariate polynomials at quadrature points.

[FEMQuad.jl](https://github.com/JuliaFEM/FEMQuad.jl): numerical integration schemes for cartesian and tetrahedron domains. I use it for obtaining Gauss-Legendre integral points in triangles, i.e. Dunavant points.

## Code Optimizations

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

Furthermore, for a system of equations with multiple variables, we can make evaluation of state variables in element `k` and quadrature point `q` into a simple linear algebra
```
ρ, ρux, ρuy, ρE = phi[:,q]' * coef_PVE[:,:,k]
```
which is equivalent to evaluate each quantity separately, but 2x faster.

Usually you don't need to worry too much about the temporary scalar allocation within loops; however, it is a good practice to check if you are not sure.

## More Things to Try

* The code I have now is based on unstructured triangular mesh. It should be pretty easy to convert that to a structured rectangular mesh, with a different set of basis functions described in the note.

* Curved boundary: not so interesting me at the moment because I am not doing high fidelity boundary layer simulations, but may be important if I have a chance to work on that in the future.

* Diffusion: DG has been practically proved to work well on advection problems, while CG is known to work well on diffusion problems. For fluid simulations, we are more interested in advection than diffusion, but we still need to include the diffusion terms. I need to be more familiar in this treatment.