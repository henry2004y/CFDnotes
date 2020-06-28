# This file was generated, do not modify it. # hide
'Computes values of pth Lagrange basis function at x between [-1,1].'
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