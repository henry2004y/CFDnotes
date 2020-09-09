@def title = "Incompressible Fluid Solver"
@def hascode = true
@def date = Date(2020, 9, 8)
@def rss = "Fluid Simulation for Computer Graphics."
@def tags = ["syntax", "code"]

\toc

# Incompressible Fluid Solver

Fluid simulation in industry focuses mostly on incompressible solvers.
Mathematically incompressible fluids have simpler expressions and can be solved numerically in an efficient way.
Here I am trying to learn some of the basic ideas and tricks, following the book Fluid Simulation for Computer Graphics, 2nd Edition and the excellent [implementation](https://github.com/henry2004y/incremental-fluids).
This is also a great opportunity for me to practice C++.

Eulerian approach is much easier in numerically approximating the spatial derivatives.

Semi-Lagrangian method is generally used for advection in solving the incompressible fluids numerically, because it's fast and often unconditionally stable.

## Euler Equation

\begin{align}
  \frac{D\mathbf{u}}{Dt} + \frac{1}{\rho}\nabla p &= \mathbf{g}\\
  \nabla\cdot\mathbf{u} &= 0
\end{align}


## Splitting Method

With operator splitting, we'll work out methods for solving these simpler equations:

\begin{align}
  \frac{Dq}{Dt} &= 0 \qquad (\text{advection}), \\
  \frac{\partial \mathbf{u}}{\partial t} &= \mathbf{g} \qquad (\text{body forces}), \\
  \frac{\partial \mathbf{u}}{\partial t} + \frac{1}{\rho}\nabla p &= 0 \qquad (\text{pressure}) \\
  \text{such that}\, \nabla\cdot\mathbf{u} &= 0 \qquad (\text{incompressibility})
\end{align}

Basic fluid algorithm:
* Start with an initial divergence-free velocity field $\mathbf{u}^0$
* For time step $n=0,1,2,...$
  * Determine a good time step $\Delta t$ to go from time $t_n$ to $t_{n+1}$.
  * Set $\mathbf{u}^A = \text{advect}(\mathbf{u}^n,\Delta t, \mathbf{u}^n)$.
  * Add $\mathbf{u}^B = \mathbf{u}^A + \Delta \mathbf{g}$.
  * Set $\mathbf{u}^{n+1} = \text{project}(\mathbf{u}^B,\Delta t)$.

## MAC Grid

In the early days of computational fluid dynamics, Harlow and Welch introduced the marker-and-cell (MAC) method method for solving incompressible flow problems.
One of the fundamental innovations of this paper was a new grid structure that makes for avery effective algorithm for enforcing incompressibility, though it may seem inconvenient for everything else.

The so-called **MAC** grid is a staggered grid, i.e., a grid where the different variables are stored at different locations. (I remember there is another name of staggered grid after people doing MHD for EM field, but it's basically the same thing.)

\fig{/assets/MAC_grid_2d.png}


## Level Set Method

This part needs more complete info!

Some geometric problems are extremely common in fluid solvers:
* When is a point inside a solid?
* What is the closest point on the surface of some geometry?
* How do we extrapolate values from one region into another?

When I wrote my simple field tracer in Julia, I spent some time on the first query especially for 1D and 2D.

The first query above suggests the right approach is **implicit surfaces**.
While there are potentially a lot of ways to generalize this, we will focus on the case where we have a continuous scalar function φ(x,y,z) which defines the geometry as follows:
* a point x is outside the geometry when φ(x) > 0,
* it’s inside the geometry when φ(x) < 0, and
* a point x is on the surface when φ(x) = 0.

In some other contexts within computer graphics, the convention may be reversed (the “outside” is the negative region) or a threshold value otherthan zero could be used, but this is usually the way we do it in fluid simulation.

**Signed distance function (SDF)**: the gradient of the function is exactly unit-length normal on the surface.

level set: a signed distance function that has been sampled on a grid.

\fig{/assets/voxel_fluid_2d.png}

A two-dimensional example of a voxelized fluid domain. **F** stands for fluid, **S** for solid, and **E** for empty.


## Advection

In the step of advection,
$$ \frac{D q}{Dt} = 0, $$
we will encapsulate this in a numerical scheme
$$ \mathbf{u}^A = \text{advect}(\mathbf{u}^n,\Delta t, q^n) $$

It bears repeating here: advection should only be called with a divergence-free velocity ﬁeld $\mathbf{u}$, i.e., one meeting the incompressibility constraint, which also satisfies the required boundary conditions.

## Projection

The heart of this scheme is to make this fluid incompressible and simultaneously enforce boundary conditions.
The implementation of
$$ \mathbf{u}^{n+1} = \text{project}(\mathbf{u}^B,\Delta t) $$
is the central topic.

The project routine will subtract the pressure gradient from the intermediate velocity field $\mathbf{u}$,
$$ \mathbf{u}^{n+1} = \mathbf{u}^n - \Delta t \frac{1}{\rho}\nabla p, $$
so that the result satisfies incompressibility inside the fluid
$$ \nabla\cdot\mathbf{u} = 0, $$
and satisfies the solid wall boundary conditions
$$ \mathbf{u}^{n+1}\cdot\widehat{n} = \mathbf{u}_{\text{solid}}\cdot\widehat{n}\qquad \text{at solid boundaries,} $$
while also respecting the free surface boundary condition that pressure be zero there:
$$ p = 0\qquad \text{at free surfaces.} $$
