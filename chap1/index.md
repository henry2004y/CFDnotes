@def title = "Math"
@def hascode = true
@def date = Date(2020, 6, 16)
@def rss = "Introduction to the required math."
@def tags = ["syntax", "code"]

\toc

# Linear Systems

## PLU Factorization

```julia:./1_2_1_PLU.jl
using LinearAlgebra
A = [4 5 8 9; 6 9 21 19; 2 1 3 4; 4 11 16 12]
F = lu(A)
@show F
```

The output is

\output{./1_2_1_PLU.jl}

## Iterative Methods

One of the most useful algorithm for solving symmetric positive semi-definite linear systems in general is called the **conjugate gradient** algorithm, usually abbreviated as CG. It's an iterative method, meaning that we start with a guess at the solution and in each iteration improve on it, stopping when we think we are accurate enough. CG chooses the iterative updates to the guess to minimize a particular measure of the error and can be guaranteed to converge to the solution eventually. Another very nice aspect of CG, as compared to Gaussian elimination for example, is that each iteration only involves multiplying A by a vector, adding vectors, multiplying vectors by scalar numbers, and computing a few dot-products—all of which are very easy to code, even in parallel, and can achieve very high efficiency on modern hardware.

The problem with CG for us, however, is that the larger the grid, the longer it takes to converge. It can be shown that the number of iterations it takes to converge to some desired accuracy is roughly proportional to the width of the grid: the maximum number of grid cells in any one direction.
In practice, when limited to a small maximum number of iterations there are simpler algorithms such as **Gauss-Seidel** and **Successive Over-Relaxation (SOR)** that tend to be much more effective than plain CG, even if they are slower at achieving full accuracy. However, there is a trick up our sleeve that can speed CG up immensely, called **preconditioning**.

**Preconditioned conjugate gradient (PCG)** is generally the desired method for symmetric positive semi-definite linear systems.
Roughly speaking, CG takes more iterations the farther $A$ is from being the identity matrix, $I$. It should be immediately obvious that solving a system with the identity matrix is pretty easy—the solution of $Ip = b$ is just $p = b$! How exactly we measure how far $A$ is from the identity involves the nature of the distribution of eigenvalues of A—but we don't need to get too detailed. The idea behind preconditioning is that the solution of $Ap = b$ is the same as the solution of $MAp = M b$ for any invertible matrix $M$. If $M$ is approximately the inverse of $A$, so that $MA$ is really close to being the identity matrix, then CG should be able to solve the preconditioned equations $MAp = M b$ really fast. PCG is just a clever way of applying CG to these preconditioned equations without actually having to explicitly form them.

### Stopping criteria

When do we know to stop? How do we check to see that our current guess is close enough to the solution? Ideally we would just measure the norm of the difference between our current guess and the exact solution—but of course that requires knowing the exact solution! So we will instead look at a vector called the **residual**:
$$ r_i = b - Ap_i.$$

That is, if $p_i$ is the ith guess at the true solution, the ith residual $r_i$ is just how far away it is from satisfying the equation $Ap = b$. When we hit the exact solution, the residual is exactly zero. Therefore, we stop our iteration when the norm of the residual is small enough, below some tolerance.

That brings us to the next question: how small is small enough? And what norm do you use to measure $r_i$? Think back to what $r_i$ means physically, for example, in the pressure solver for the incompressible fluid. These equations resulted from deriving that $b − Ap$ is the negative of the finite difference estimate of the divergence of $\mathbf{u}^{n+1}$, which we want to be zero. Thus the residual measures exactly how much divergence there will be in the velocity field after we've updated it with our current estimate of the pressure. It seems sensible, then, to take the infinity norm of the residual (the maximum absolute value of any entry) and compare that to some small number $tol$, so that we know the worst case of how compressible our new velocity field could possibly be. Often in practice we just pick an arbitrary small fixed number for $tol$, say $10^{-6}\,\ text{s}^{-1}$.
Smaller tolerances will result in less erroneous divergence but will take more iterations (and time) to compute, so there's a clear trade-off in adjusting this number up or down. This is particularly important since typically most of the time in these kind of fluid simulations is spent in the pressure solver: this is the code that usually demands the most optimizing and tuning.

That said, the fact that our absolute tolerance on stopping PCG has physical units (of one over time) can also be worrisome.
It may also be desired to have a generic linear solver available in the code, which doesn’t make assumptions about the units of the system.
Therefore we often use a relative residual measure for convergence: stop when the ratio between the current residual norm and the initial right-hand-side norm is smaller than a given dimensionless tolerance (with no physical units needed). That is, stop when
$$ ||r|| = ||b-Ap|| \le tol||b||. $$ 

Furthermore, we have to deal with the inexact floating point calculations and the computation time constraint.

### Initial guess

One nice thing about PCG is that if we start with a good guess, we can get to an acceptable solution much faster.
There are generally two options:
1. using the solution from last time,
2. all zeros.

In the more interesting/dynamic cases (and really, why would you be simulating fluids that are just sitting still?), the pressure can change significantly from one time step to the next, or in fact may be defined on different grid cells anyhow (e.g., as liquid moves from one grid cell to another) and can’t be used. Therefore we usually use the vector of all zeros as the initial guess.

### BLAS

Using highly optimized libraries is recommended. **Basic Linear Algebra Subroutine (BLAS)** lies underneath literally all of the modern programming languages for scientific computing.

### Preconditioner

From the stand-point of convergence the perfect preconditioner would be $A^{−1}$, except that’s obviously far too expensive to compute. The true ideal is something that is both fast to compute and apply, and is effective in speeding up convergence, so as to minimize the total solution time.

#### Incomplete Cholesky

There are many, many choices of preconditioner, with more being invented each year. 
One classical preconditioner used in solving the pressure equation in incompressible fluid is from the **Incomplete Cholesky** (IC) family, which is highly robust in handling irregular domains (like the shape of a liquid splash). Its chief problems are that it’s hard to parallelize effectively and that it's not optimally scalable (the number of iterations required for PCG slowly increases with grid size. Algorithms in the **domain decomposition** and **multigrid** family of methods can provide both excellent parallelization and optimal scalability (the time it takes to solve a problem is linearly proportional to thenumber of grid cells) but are not trivial to implement in a way that is robust to irregular domains.

Here we will discuss a simple form of domain decomposition that can provide good parallelization for Incomplete Cholesky.
Recall how you might directly solve a system of linear equations with Gaussian elimination. That is, you perform row reductions on the system until the matrix is upper-triangular, and then use back substitution to get the solution one entry at a time. Mathematically, it turns out that this is equivalent to factoring the matrix A as the product of a lower- and an upper-triangular matrix and then solving the triangular systems one after the other.
In the case of a symmetric positive definite A, we can actually do it so that the two triangular matrices are transposes of each other:
$$ A = LL^T. $$
This is called the Cholesky factorization. The original system $Ap = b$ is the same as $L(L^Tp) = b$, which we can solve as
\begin{align}
  &\text{solve}\quad Lq = b\quad &\text{with forward substitution},\\
  &\text{solve}\quad L^Tp = q\quad &\text{with backward substitution}.
\end{align}
The main reason that we don't typically use this method for fluids is that although A has very few non-zeros, L can have a lot.
In three dimensions the amount of **fill-in** (extra non-zeros in L) is particularly bad; **direct solvers** that take this approach can easily fail on 3D problems due to lack of memory.

Basic Incomplete Cholesky tackles this problem with a very simple idea: whenever the Cholesky algorithm tries to create a new non-zero in a location that is zero in A, cancel it—keep it as zero. On the one hand, the resulting L is just as sparse as A, and memory is no longer an issue. On the other hand, we deliberately made errors: $A\neq LL^T$ now. However, hopefully the "incomplete" factorization is close enough to A that doing the solves in Eq.(4) is close enough to applying $A^{−1}$ so that we have a useful preconditioner for PCG!

Technically, performing Incomplete Cholesky only allowing non-zeros in $L$ where there are non-zeros in $A$ is called level zero: IC(0). There are variations that allow a limited number of non-zeros in other locations. For the relatively simple Laplacian matrix we are dealing with, they generally are not worth the computational effort.

To make this more precise, IC(0) constructs a lower-triangular matrix $L$ with the same non-zero pattern as the lower triangle of $A$, such that $LL^T= A$ in the locations where $A$ is non-zero. The only error is that $LL^T$ is non-zero in some other locations where A is zero.

Assume we order our grid cells (and the corresponding rows and columns of $A$) lexicographically, say along the i-dimension first, then the j-dimension, and finally the k-dimension. Suppose we split $A$ up into its strict lower triangle $F$ and diagonal $D$:
$$ A = F + D + F^T. $$
Then, it can be shown for the particular $A$ we're solving, IC(0) factor $L$ is of the form
$$ L = F E^{−1}+ E, $$
where $E$ is a diagonal matrix. That is, all we need to compute and store are the diagonal entries of $L$, and we can infer the others just from $A$!

Crunching through the algebra gives the following formulas for computing the entries along the diagonal of $E$. In two dimensions,
$$ E_{(i,j)}=\sqrt{A_{(i,j),(i,j)} − (A_{(i−1,j),(i,j)}/E_{(i−1,j)})^2 − (A_{(i,j−1),(i,j)}/E_{(i,j−1)})^2}. $$
In three dimensions,
$$ E_{(i,j,k)}=\sqrt{A_{(i,j,k),(i,j,k)} − (A_{(i−1,j,k),(i,j,k)}/E_{(i−1,j,k)})^2 − (A_{(i,j−1,k),(i,j,k)}/E_{(i,j−1,k)})^2 − (A_{(i,j,k−1),(i,j,k)}/E_{(i,j,k−1)})^2} $$

In these equations, we replace terms referring to a non-fluid cell (or cell that lies off the grid) with zero. 

#### Modified Incomplete Cholesky

Incomplete Cholesky is a good preconditioner that can effectively reduce our iteration count when solving the pressure equations, and is often the default choice when preconditioning any general matrix. 
But, for almost no extra cost, we can do better for our particular $A$! A slight tweak to IC, *Modified* Incomplete Cholesky (MIC), scales significantly better: if our grid is $n$ grid cells wide, regular IC(0) will require $O(n)$ iterations but MIC(0) will converge in only $O(n^{1/2})$ iterations, with a fairly low hidden constant. Modified Incomplete Cholesky works exactly like Incomplete Cholesky, except instead of just discarding those unwanted non-zeros, we account for them by adding them to the diagonal of $L$.

To make this more precise, MIC(0) constructs a lower-triangular matrix $L$ with the same non-zero pattern as the lower triangle of $A$, such that
* The off-diagonal non-zero entries of $A$ are equal to the corresponding ones of ($LL^T$).
* The sum of each row of $A$ is equal to the sum of each row of ($LL^T$).

This boils down to a slightly different calculation for the diagonal entries: the modified $L$ is also equal to $F E^{−1} + E$, just for a different $E$. In two dimensions,
$$ E_{(i,j)} = \sqrt{ A_{(i,j),(i,j)} − (A_{(i−1,j),(i,j)}/E_{(i−1,j)})^2 − (A_{(i,j−1),(i,j)}/E_{(i,j−1)})^2 − A_{(i−1,j),(i,j)}A_{(i−1,j),(i−1,j+1)}/E^2_{(i−1,j}) − A_{(i,j−1),(i,j)}A_{(i,j−1),(i+1,j−1)}/E^2_{(i,j−1)} }. $$
In three dimensions,
$$ E(i,j,k)=\sqrt{ A_{(i,j,k),(i,j,k)} − (A_{(i−1,j,k),(i,j,k)}/E_{(i−1,j,k)})^2 − (A_{(i,j−1,k),(i,j,k)}/E_{(i,j−1,k)})^2 − (A_{(i,j,k−1),(i,j,k)}/E_{(i,j,k−1)})^2 − A_{(i−1,j,k),(i,j,k)} (A_{(i−1,j,k),(i−1,j+1,k)} + A_{(i−1,j,k),(i−1,j,k+1)})/E^2_{(i−1,j,k)} − A_{(i,j−1,k),(i,j,k)} (A_{(i,j−1,k),(i+1,j−1,k)} + A_{(i,j−1,k),(i,j−1,k+1)})/E^2_{(i,j−1,k)} − A_{(i,j,k−1),(i,j,k)} (A_{(i,j,k−1),(i+1,j,k−1)} + A_{(i,j,k−1),(i,j+1,k−1)})/E^2_{(i,j,k−1)} } $$

If you're curious, the intuition behind MIC (and why it outperforms IC) lies in a Fourier analysis of the problem. If you decompose the error as a superposition of Fourier modes, some low frequency (smooth) and some high frequency (sharp), it turns out IC is only effective at removing the high-frequency components of error. On the other hand, MIC is forced to match the action of $A$ on the lowest frequency mode of all, the constant, and thus is more effective at all frequencies. [^1]

[^1]: Continuing this train of thought, looking for methods that work well on all frequency components of the error can lead to multigrid that explicitly solves the equations at multiple resolutions.

In practice, you can squeeze out even better performance by taking a weighted average between the regular Incomplete Cholesky formula andthe modified one, typically weighting with 0.97 or more (getting closer to 1 for larger grids).

Additionally we can have a built-in safety tolerance. In some situations, such as a single-cell–wide line of fluid cells surrounded by solids, IC(0) and MIC(0) become exact—except that $A$ is singular in this case: the exact Cholesky factorization doesn't exist. This is manifested by hitting a zero or—when rounding error is factored in—very small value for $e$, and it can safely be cured by replacing that small value with, for example, the diagonal entry from $A$. This safety check also comes in handy if you want to solve more general linear systems, where the Incomplete Cholesky factorization (and even more so the Modified Incomplete Cholesky factorization) may fail to exist without this check.

#### Domain decomposition

Solving the linear system for pressure is one of the more expensive steps in a typical fluid solve, so parallelization is crucial. Unfortunately, the computation of Incomplete Cholesky described above, running through the system in lexicographic order, is inherently sequential: the forward and backward triangular solves needed at every step of PCG have data dependency that eliminates the possibility for parallelization.
It is possible to reorder thematrix and the grid cells in a way that allows a parallel factorization, but interestingly and frustratingly enough, such an ordering seriously harms the power of the preconditioner, rendering it not much better than doing nothing at all. Research has continued in this vein using *graph multicolorings* for parallelism, but generally a more complex version of Incomplete Cholesky with more nonzeros in the sparsity pattern is required at a minimum.

Luckily, we can take a different approach to parallelization, domain decomposition. Domain decomposition is essentially the divide-and-conquer principle applied to solving sparse linear systems: we partition our "domain" (the grid of unknowns) into "subdomains" (subsets of the unknowns), solve the smaller linear system restricted to each subdomain independently, and patch the results together into an approximate solution for the original domain. Ideally almost all of the work takes place in solving the subdomains, which can all be done in parallel, leading to excellent parallel efficiency.

Here we will go with a particularly simple variant of what is termed "Additive Overlapping Schwarz".

First we need to partition our grid into subdomains. For simplicity, we may as well take them to be subgrids. The critical requirements are:
* We need at least as many subgrids as parallel threads we want to use (fewer is better for convergence, but having more than one subgrid per thread can make dynamic load balancing work better).
* The subgrids, taken all together, have to cover the entire original grid: we can't leave any gaps.
* The subgrids should overlap a bit, i.e. have some grid cells in common with their neighbors. The bigger the overlap, the fewer the number of iterations we will need but the slower to compute each iteration will be (and the more memory required).

For example, if we had a $100\times 100\times 100$ original grid, we could split it into eight $52\times 52\times 52$ subgrids, with an overlap band of four grid cells—the grid cells in the very center of the domain would be a part of all the subgrids.

For each call to the preconditioner, we have to first solve or approximately solve each subdomain's part of the equation $Az = r$ independently. For each subdomain, this entails taking just the entries of $A$ and $r$ in the grid cells of that subdomain, and ignoring the offdiagonal entries of $A$ that reference grid cells outside the subdomain. The solve could be a full linear solve to maximum precision using some other technique, but it is likely to be more efficient to just use one call to MIC(0) instead (giving a very approximate solution). Because this is going to part of an outer PCG loop, not every method is allowed here—approximations that are unsymmetric, or nonlinear (such as taking a few "inner" iterations of PCG without converging fully) can cause convergence troubles. 
Because the solves are completely independent, this step is trivial to parallelize without synchronization or communication required: just some shared reads to $A$ and $r$ (but *not* shared writes). Once each subdomain has its own local estimate of the solution $z$ restricted to its subdomain, we simply add up all the subdomain solves into one global estimate $z$, which is returned to the outer PCG loop. This addition, in the overlap regions where multiple subdomains contribute, does require some form of synchronization or locking—or the calculation of the overlap entries can be left to a single thread at the end.

When the number of subdomains becomes large, the convergence rate of plain domain decomposition diminishes. To get an optimal high-performance preconditioner, it is necessary to include a **coarse grid correction**. This is basically a step (similar to multigrid) inside the preconditioner where we solve a very coarse version of the problem, with one or just a few variables per subdomain, representing the smoothest components of the error. [^2]

[^2]: The book Discretely-Discontinuous Galerkin (DDG) Method of Edwards and Bridson provides a very simple implementation that can be applied to just about any sort of PDE problem, while providing extremely good speed-up.

Packages exist, like [LimitedLDLFactorizations.jl](https://github.com/JuliaSmoothOptimizers/LimitedLDLFactorizations.jl)

## Krylov Subspace Methods

There is a package [IterativeSolvers.jl](https://juliamath.github.io/IterativeSolvers.jl/dev/) that provides efficient iterative algorithms for solving large linear systems, eigenproblems, and singular value problems.
Most of the methods can be used matrix-free.

For example,

```julia:./1_2_3_GMRES.jl
using IterativeSolvers
A = rand(3,3)
A = A'*A
b = rand(3)
x = gmres(A,b)
println(x)
```

gives

\output{./1_2_3_GMRES.jl}

# Numerical Integration

## Basic polynomials

[Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) provides basic arithmetic, integration, differentiation, evaluation, and root finding over dense univariate polynomials.

[SpecialPolynomials.jl](https://github.com/jverzani/SpecialPolynomials.jl) provides many more classic orthogonal polynomials, the Lagrange and Newton interpolating polynomials, and the Bernstein polynomials.

## Legendre polynomials

There are packages like [PolynomialBases.jl](https://github.com/ranocha/PolynomialBases.jl).

For example, to evaluate the `p` order Legendre polynomial at `x`,

```julia:./1_3_legendre.jl
using PolynomialBases
x = 0.0
p = 1
@show legendre(x,p)
```

outputs

\output{./1_3_legendre.jl}

## Lagrange polynomials

To evaluate the value of Lagrange basis function on [-1,1]:

```julia:./lagrange.jl
"Computes values of pth Lagrange basis function at x between [-1,1]."
function lagrange(p::Integer, x)

   @assert p ≥ 0

   n = p + 1
   ϕ = ones(eltype(x),length(x),n)
   if n == 1
      return phi
   end

   InterPts = range(-1, 1, length=n)

   for i = 1:n
      for j = 1:n
         if j != i
            ϕ[i] *= (x - InterPts[j]) / (InterPts[i] - InterPts[j])
         end
      end
   end
   return ϕ
end

@show lagrange(1, 0.0)
```

gives
\output{./lagrange.jl}


A even easier approach would be taking advantage of the [Polynomials.jl](https://github.com/JuliaMath/Polynomials.jl) package and return the coefficients of the polynomials:

```julia:./lagrangeBasis.jl
using Polynomials

function lagrangeBasis(p::Integer)
   @assert p ≥ 0

   p == 0 && return [one(Polynomial)]

   roots = range(1, -1, length=p+1)
   ϕ = fill(one(Polynomial),p+1)
   ψ = Vector{Polynomial}(undef,p+1)
   for i = 1:p+1
      for j = 1:p+1
         if i != j
            temp = -1.0 / (roots[i]-roots[j])
            ψ[j] = Polynomial([roots[j]*temp, temp])
            ϕ[i] *= ψ[j]
         end
      end
   end
   return ϕ
end

@show lagrangeBasis(1)
@show lagrangeBasis(1)[1](1.0) # evaluate the first basis at 1.0
```

outputs

\output{./lagrangeBasis.jl}


## Barycentric polynomials

Another version of Lagrange polynomials, but more stable.
One implementation can be found at [BarycentricInterpolation.jl](https://github.com/dawbarton/BarycentricInterpolation.jl)

# Time Marching

Many functionalities, especially for ODEs, can be found at [DifferentialEquations.jl](https://docs.sciml.ai/stable/).

## CFL Condition

Each point of a numerical solution also has a domain of dependence: again, the set of locations in the initial conditions that have an effect on the value of the solution at that point.
It should be intuitively obvious that the numerical domain of dependence, at least in the limit, must contain the true domain of dependence if we want to get the correct answer. This is, in fact, the CFL condition: convergence is only possible in general if, in the limit as $\Delta x \rightarrow 0$ and $\Delta t \rightarrow 0$, the numerical domain of dependence for each point contains the true domain of dependence.

## Splines

中文对应样条插值。There are package like [Dierckx.jl](https://github.com/kbarbary/Dierckx.jl), [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl), and [GridInterpolations.jl](https://github.com/sisl/GridInterpolations.jl) for this purpose.

## Discretization Approaches
