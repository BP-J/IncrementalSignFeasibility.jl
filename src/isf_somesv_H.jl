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
`isf_somesv_H!(V, perm, info, options)`

Modifies info during the computation of (very) few stem vectors.
Compute the stem vectors of V that can be deduced from the r linear 
independent vectors given by the QR factorization and put them in the 
matrix info.stems.
"""
function isf_somesv_H!(V::Matrix, perm::Vector{Int}, info::Info, options::Options)
    p = size(V,2)
    r = info.r
    info.nb_stems_sym = p-r                                         # info.nb_stems_sym equals p-r by definition
    
    ## computation of the p-r stem vectors given by the QR factorization: 
    # the first r vectors are supposed independent (after use of the permutation)
    # so the remaining p-r can be expressed as combinations of the first
    # these naturally form some stem vectors

    for ivp in perm[r+1:p]
        cols = [perm[1:r];ivp]                                      # considered columns
        Z = nullspace(V[:,cols])                                    # element in the null space
        I = findall(x -> abs(x) > options.tol_coordinates, Z)       # detection of nonzero elements
        stem = zeros(Int,p)
        stem[cols[I]] = sign(Z[I[1]])*sign.(Z[I])                   # getting the signs
        append!(info.stems_sym_init, stem)                          
    end

    # the stem vectors are necessarily different, so no 'unique'-type call needed

    info.stems_sym = info.stems_sym_init                            # the Elastic part will disappear by itself
    info.stem_sizes_sym = vec(sum(abs.(info.stems_sym), dims=1))    # the sum returns a 1 x [size] matrix, not a vector
    # nothing to do for the asym ones, there are none

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
        info.stem_zero_indices_sym[1] = findall(iszero, info.stems_sym[perm[p], 1:info.nb_stems_sym] )

        for j in 2:(p-r-1)
            I = info.stem_zero_indices_sym[j-1]
            info.stem_zero_indices_sym[j] = I[ findall(iszero, info.stems_sym[perm[p+1-j], I]) ]
        end
    
    end
end