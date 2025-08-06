#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Version 0.1, 2025.
#
# If you found this piece of software useful, please cite the paper:
# J.-P. Dussault, J.Ch. Gilbert, B. Plaquevent-Jourdain,
# "Primal and Dual Approaches for the Chamber Enumeration
# of real hyperplane arrangements", 2025.
#
# Authors:
# - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
# - Jean Charles Gilbert (INRIA, France),
# - Baptiste Plaquevent-Jourdain (INRIA & Univ. of Sherbrooke, Canada).
#
# Copyright 2025, INRIA (France) and Univ. of Sherbrooke (Canada).
#
# ISF is distributed under the terms of the Q Public License version
# 1.0.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Q Public
# License version 1.0 for more details.
#
# You should have received a copy of the Q Public License version 1.0
# along with this program. If not, see
# <https://doc.qt.io/archives/3.3/license.html>.
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

"""
`isf_somesv_HnH!(V, perm, info, options)`

Modifies info during the computation of (very) few stem vectors.
Compute the stem vectors of V that can be deduced from the r linear 
independent vectors given by the QR factorization and put them in the 
matrix info.stems.
"""
function isf_somesv_HnH!(Vt::Matrix, perm::Vector{Int}, info::Info, options::Options)
    p = size(Vt,2)
    n = size(Vt, 1) - 1
    r = info.r
    
    ## computation of the p-r stem vectors given by the QR factorization: 
    # the first r vectors are supposed independent (after use of the permutation)
    # so the remaining p-r can be expressed as combinations of the first
    # these are naturally some stem vectors

    for ivp in perm[r+1:p]
        cols = [perm[1:r];ivp]                                      # considered columns
        Z = nullspace(Vt[1:n,cols])                                 # element in the null space
        I = findall(x -> abs(x) > options.tol_coordinates, Z)       # indices of the nonzero components - the real circuit
        stem = zeros(Int, p)
        if abs(Vt[n+1,cols]'*Z[:]) <= options.tol_nonzero_q         # considered in Vt's nullspace, so "asymmetric" stem vector
            stem[cols[I]] = sign(Z[I[1]])*sign.(Z[I])               # one could also use the sign(T*...) convention but not necessary
            append!(info.stems_asym_init, stem)
            info.nb_stems_asym += 1
        else                                                        # so "symmetric" stem vector
            stem[cols[I]] = sign(Vt[n+1,cols]'*Z[:])*sign.(Z[I])    # sign convention, see above
            # the scalar product cannot be close to 0, otherwise Z would be in 
            # the nullspace of Vt which is the previous case.
            append!(info.stems_sym_init, stem)
            info.nb_stems_sym += 1                    
        end
    end

    # the stem vectors are necessarily different, so no 'unique'-type call needed

    info.stems_sym = info.stems_sym_init                            # the Elastic part will disappear by itself
    info.stems_asym = info.stems_asym_init                          # the Elastic part will disappear by itself
    info.stem_sizes_sym  = vec(sum(abs.(info.stems_sym), dims=1))   # the sum returns a 1 x [size] matrix, not a vector
    info.stem_sizes_asym = vec(sum(abs.(info.stems_asym), dims=1))  # the sum returns a 1 x [size] matrix, not a vector
    
    ### in case the heuristic C is not used, a particular reordering of the stem vectors is used.
    # Said quickly, the idea is that, in this case, the order of the vectors added is fixed by perm.
    # Thus, the indices will be added in the order chosen by perm (thus through the QR factorization).
    # Then, consider a stem vector is given (the permutation forgotten so simplify the notation) by
    # 1-2-3; then one would like to compare a given sign vector s^k = (s_1, ... s_k) with it before
    # comparing with a stem vector that has indices p-2 p-1 p for instance.

    # The first array in the list corresponds to the ones having a 0 in last (taking perm into account)
    # position, the second array the ones with zeros in last two positions and so on

    if options.bestv == 0 && p > r+1
        info.stem_zero_indices_sym = [Any[] for _ in 1:(p-r-1)]     # the first 2 ones are kept empty, as stem vectors have sizes > 2 
        info.stem_zero_indices_asym = [Any[] for _ in 1:(p-r)]      # the first   one   is kept empty, as stem vectors have sizes > 1

        info.stem_zero_indices_sym[1] = findall(iszero, info.stems_sym[perm[p], 1:info.nb_stems_sym] )
        info.stem_zero_indices_asym[1] = findall(iszero, info.stems_asym[perm[p], 1:info.nb_stems_asym] )

        for j in 2:(p-r-1)
            I = info.stem_zero_indices_sym[j-1]
            info.stem_zero_indices_sym[j] = I[ findall(iszero, info.stems_sym[perm[p+1-j], I]) ]
            I = info.stem_zero_indices_asym[j-1]
            info.stem_zero_indices_asym[j] = I[ findall(iszero, info.stems_asym[perm[p+1-j], I]) ]
        end
    
    end
end