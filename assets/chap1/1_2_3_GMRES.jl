# This file was generated, do not modify it. # hide
using IterativeSolvers
A = rand(3,3)
A = A'*A
b = rand(3)
x = gmres(A,b)
println(x)