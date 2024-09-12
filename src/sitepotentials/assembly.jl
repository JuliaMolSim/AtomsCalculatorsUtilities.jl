

# ---------------------------------------------------
# utilities 

function site_virial(dV, Rs) 
   TV = promote_type(eltype(eltype(dV)), eltype(eltype(Rs)))
   return - sum( dv_i * 𝐫_i' for (dv_i, 𝐫_i) in zip(dV, Rs); 
                 init = zero(SMatrix{3, 3, TV}) )
end



# ---------------------------------------------------
# main assembly codes 

AtomsCalculators.@generate_interface function AtomsCalculators.potential_energy(
                  sys, 
                  V::SitePotential; 
                  domain = 1:length(sys), 
                  executor = ThreadedEx(), 
                  nlist = PairList(sys, cutoff_radius(V)), 
                  kwargs...)
   uE = energy_unit(V)
   E = Folds.sum( domain, executor; 
                  init = zero_energy(sys, V) 
   ) do i
      Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, i) 
      e_i = eval_site(V, Rs, Zs, z0) * uE 
   end

   return E
end


AtomsCalculators.@generate_interface function AtomsCalculators.virial(
                  sys, 
                  V::SitePotential; 
                  domain   = 1:length(sys), 
                  executor = ThreadedEx(),
                  nlist    = PairList(sys, cutoff_radius(V)),
                  kwargs...)
   uE = energy_unit(V)
   vir = Folds.sum( domain, executor; 
                    init = zero_virial(sys, V) 
   ) do i 
      Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, i)
      Ei, ∇Ei = eval_grad_site(V, Rs, Zs, z0)
      vir_i = site_virial(∇Ei, Rs) * uE
   end

   return vir
end


function energy_forces_virial(
                  sys, 
                  V::SitePotential; 
                  domain   = 1:length(sys), 
                  executor = ThreadedEx(),
                  ntasks   = Threads.nthreads(),
                  nlist    = PairList(sys, cutoff_radius(V)),
                  kwargs...)
   uE = energy_unit(V)
   fE = force_unit(V)
   E_F_V = Folds.sum( collect(chunks(domain, ntasks)), 
                      executor;
                      init=[ zero_energy(sys, V), 
                             zero_forces(sys, V), 
                             zero_virial(sys, V) ]
   ) do (sub_domain, _)

      E = zero_energy(sys, V)
      frc = zero_forces(sys, V)
      vir = zero_virial(sys, V)

      for i in sub_domain
         Js, Rs, Zs, z0 = get_neighbours(sys, V, nlist, i)
         Ei, ∇Ei = eval_grad_site(V, Rs, Zs, z0)
         E += Ei * energy_unit(V)
         vir += site_virial(∇Ei, Rs) * energy_unit(V)
         # we have to update the forces in a loop since Js indices can be repeated
         for α in 1:length(Js)
            frc[Js[α]] -= ∇Ei[α] * force_unit(V)
          end
          frc[i] += sum(∇Ei) * force_unit(V)
      end
      [E, frc, vir] 
   end
   
   return (energy = E_F_V[1], forces = E_F_V[2], virial = E_F_V[3],)
end


function AtomsCalculators.energy_forces(at, V::SitePotential; kwargs...)
   efv = energy_forces_virial(at, V; kwargs...)
   return (energy = efv.energy, forces = efv.forces)
end

AtomsCalculators.@generate_interface function AtomsCalculators.forces(at, V::SitePotential; kwargs...) 
   return AtomsCalculators.energy_forces(at, V; kwargs...)[:forces]
end
