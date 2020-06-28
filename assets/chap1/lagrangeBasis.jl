# This file was generated, do not modify it. # hide
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