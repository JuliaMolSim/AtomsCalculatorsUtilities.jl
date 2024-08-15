
using NeighbourLists
using NeighbourLists: PairList 


function _neigs_range(nlist::PairList, i0) 
   n1, n2 = nlist.first[i0], nlist.first[i0+1]-1
   return n1:n2
end

function get_neighbours!(out, sys, V, nlist::PairList, i0)
   nrange = _neigs_range(nlist, i0)
   Js = out.Js; Rs = out.Rs; Zs = out.Zs
   for (n, j) in enumerate(nrange)
      Js[n] = nlist.j[j]
      Rs[n] = NeighbourLists._getR(nlist, j)
      Zs[n] = get_id(sys, V, Js[n])
   end
   z0 = get_id(sys, V, i0) 
   return Js, Rs, Zs, z0
end

function whatalloc(::typeof(get_neighbours!), 
                   sys, V, nlist::PairList, i0)
   num_neigs = length(_neigs_range(nlist, i0))
   TR = eltype(nlist.X)
   TZ = typeof(get_id(sys, V, 1))
   return (Js = (Int, num_neigs), 
           Rs = (TR, num_neigs),
           Zs = (TZ, num_neigs), )
end


"""
`get_neighbours(sys, V, nlist, i) -> Js, Rs, Zs, z0`
"""
function get_neighbours(sys, V, nlist::PairList, i)
   # allocate the output arrays 
   alc = whatalloc(get_neighbours!, sys, V, nlist, i)
   out = ( Js = zeros(alc.Js...), 
           Rs = zeros(alc.Rs...), 
           Zs = zeros(alc.Zs...) )
   # do the actual computation            
   return get_neighbours!(out, sys, V, nlist, i)           
end

# the old reference implementation for testing purposes. 
function get_neighbours_old(sys, V, nlist::PairList, i) 
   Js, Rs = NeighbourLists.neigs(nlist, i)
   Zs = [ get_id(sys, V, j) for j in Js ]
   z0 = get_id(sys, V, i) 
   return Js, Rs, Zs, z0 
end 