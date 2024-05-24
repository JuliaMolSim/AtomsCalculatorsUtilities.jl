
import ForwardDiff

"""
Implements `hessian_site` using ForwardDiff on the site gradient. 
This does NOT automatically implement the hessian for site potential. 
A site potential implementation  must overload `hessian_site` e.g. as follows
```julia
hessian_site(V::StillingerWeber, Rs::AbstractVector{SVector{D, TF}}, Zs, z0
                   ) where {D, TF} = ad_hessian_site(V, Rs, Zs, z0)
```
"""
function ad_hessian_site(V::SitePotential, 
                               Rs::AbstractVector{SVector{D, TF}}, 
                               Zs, z0) where {D, TF}
   Rs2x = Rs -> reinterpret(eltype(eltype(Rs)), Rs)
   x2Rs = x -> reinterpret(SVector{D, eltype(x)}, x)
   # f = gradient in x coordinates 
   f = x -> Rs2x( eval_grad_site(V, x2Rs(x), Zs, z0)[2] ) 
   # ∂f = jacobian of f so i.e. the hessian of site energy in x coords
   ∂f = ForwardDiff.jacobian(f, Rs2x(Rs))
   return ∂f
end


"""
Implements `block_hessian_site` using ForwardDiff on the site gradient. 
This does NOT automatically implemet the hessian. A site potential implementation 
must overload `block_hessian_site` e.g. as follows
```julia
block_hessian_site(V::StillingerWeber, Rs::AbstractVector{SVector{D, TF}}, Zs, z0
                   ) where {D, TF} = ad_block_hessian_site(V, Rs, Zs, z0)
```
"""
function ad_block_hessian_site(V::SitePotential, 
                               Rs::AbstractVector{SVector{D, TF}}, 
                               Zs, z0) where {D, TF}
   ∂f = ad_hessian_site(V, Rs, Zs, z0)                               
   # convert to smatrix blocks 
   nR = length(Rs)
   Hblock = zeros(SMatrix{D, D, TF}, nR, nR)
   for j1 = 1:nR, j2 = 1:nR 
      J1 = 3 * (j1-1) .+ (1:3)
      J2 = 3 * (j2-1) .+ (1:3)
      Hblock[j1, j2] = SMatrix{D, D}(∂f[J1, J2])
   end
   return Hblock 
end


# ------------------------------------------------ 
#   global hessian assembly prototypes 

function hessian(sys, V::SitePotential; 
                       domain=1:length(sys), 
                       executor=ThreadedEx(), 
                       nlist=nothing, 
                       kwargs...
                       )
   if isnothing(nlist)
      nlist = PairList(sys, cutoff_radius(V))
   end

   Nat = length(sys)
   D = n_dimensions(sys)

   # learn what the types are 
   Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, 1) 
   E1 = eval_site(V, Rs, Zs, z0) 
   TF = typeof(E1)
   hU = energy_unit(V) / length_unit(V)^2
   TFU = typeof(E1 * hU)
   # assuming here that the type of everything else will be the same?
   # need to think whether this is a restriction. 
   # also need to do something about fucking units 

   H = zeros(TFU, D*Nat, D*Nat)

   for i in domain 
      Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, i) 
      Hi = hessian_site(V, Rs, Zs, z0)
      release!(Js); release!(Rs); release!(Zs)

      Nr = length(Js)
      Ji = (i - 1) * D .+ (1:D)
      for (α1, j1) in enumerate(Js), (α2, j2) in enumerate(Js)
         A1 = (α1-1) * D .+ (1:D)
         A2 = (α2-1) * D .+ (1:D)
         J1 = (j1-1) * D .+ (1:D)
         J2 = (j2-1) * D .+ (1:D)
         H[J1, J2] += Hi[A1, A2] .* hU
         H[J1, Ji] -= Hi[A1, A2] .* hU
         H[Ji, J2] -= Hi[A1, A2] .* hU
         H[Ji, Ji] += Hi[A1, A2] .* hU
      end
   end

   return H
end





function block_hessian(sys, V::SitePotential; 
                       domain=1:length(sys), 
                       executor=ThreadedEx(), 
                       nlist=nothing, 
                       kwargs...
                       )
   if isnothing(nlist)
      nlist = PairList(sys, cutoff_radius(V))
   end

   Nat = length(sys)
   D = n_dimensions(sys)

   # learn what the types are 
   Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, 1) 
   E1 = eval_site(V, Rs, Zs, z0) 
   TF = typeof(E1)
   hU = energy_unit(V) / length_unit(V)^2
   TFU = typeof(E1 * hU)

   # assuming here that the type of everything else will be the same?
   # need to think whether this is a restriction. 
   # also need to do something about fucking units 

   H = zeros(SMatrix{D, D, TFU}, Nat, Nat)

   for i in domain 
      Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, i) 
      Hi = block_hessian_site(V, Rs, Zs, z0)
      release!(Js); release!(Rs); release!(Zs)

      # Hi[j1, j2] = ∂²V/∂𝐫_ij1 ∂𝐫_ij2
      # where 𝐫_ij1 = 𝐫_j1 - 𝐫_i
      nRs = length(Js)
      for (α1, j1) in enumerate(Js), (α2, j2) in enumerate(Js)
         H[j1, j2] += Hi[α1, α2] * hU
         H[j1, i] -= Hi[α1, α2] * hU
         H[i, j2] -= Hi[α1, α2] * hU
         H[i, i] += Hi[α1, α2] * hU
      end
   end

   return H
end

