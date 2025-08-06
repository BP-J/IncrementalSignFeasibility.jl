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
`isf_allsv_HnH!(Vt, perm, info, options, values)`

Modifies info by calling `isf_allsv_from_indices_HnH!` on subsets of [1:p].
Computes all stem vectors of V | Vt and puts them in info.
It proceeds by looking at all the potential subsets of columns, then
computes their nullity; if it is 0, the vectors are independent
if it is 1, then the vector of the null space (up to constant)
gives the corresponding stem vector(s). 

The method is a bit more tedious due to the nature of the computations. 
In particular, if a subset finds a circuit for V, then it may be one for 
Vt as well (unlikely), otherwise a circuit containing this one is searched
in R^{n+1}: the dimension is n + HnH_info, which tells the status 
of the recurrence in the local subbranch. 

Note that the maximal subset size is r+2 [= (r+1) + 1 for the n+1th dimension]
instead of n+1 previously and the minimal size 2. 

WARNING: the "asymmetric" stem vectors here correspond to the stem vectors 
of Vt in dimension n+1. These are the "symmetric" stem vectors in the nH situation. 

More precisions are found in `isf_allsv_from_indices_HnH.jl`.
The 'symmetric' (H) and 'asymmetric' (nH) versions are in their respective 
functions (files) `isf_allsv_H!(.jl)`, `isf_allsv_nH!(.jl.)`.

The computations can be done differently (and often slightly faster), see
`isf_allsv_HnH_Wechelon(_rec).jl`.
"""
function isf_allsv_HnH!(Vt::Matrix, perm::Vector{Int64}, info::Info, options::Options, values::Values)

    p = size(Vt,2)
    if p < 3
        println("isf_allsv_HnH!: function is not made for small values of p = $(p)")
        info.flag = values.fail_technicality
    end
    r = info.r

# ----------------------------------------------------------------------------------------------------
# Compute all the stem vector of a "matrix" without permuting the columns (version with a binary tree)
# ----------------------------------------------------------------------------------------------------
    
    # A difference is considering the subsets of size 2, since two parallel hyperplanes
    # correspond to a circuit but the hyperplanes are not identical.

    # since every subset of 2 <= size <= r+1 is considered, the first one is the one with vectors 1,2
    bvec = zeros(Int, p)
    bvec[1:2] = [1,1]                                                                           # this one
    p_bvec = 2                                                                                  # binary tree pointer
    HnH_info = 0                                                                                # info: dimension = n + HnH_info
    memory = 1                                                                                  # info: backtracking index, to know where/when the recursion stays in dimension n+1
    type = 0                                                                                    # info: stem type (defined here to avoid definition error)

    if options.rational                                                                         # exact (rational) computations; test before to avoid doing it every time
        while p_bvec > 0                                                                        # the loop ends when the pointer is zero
            newstem = false
            if 1 < sum(bvec) < r+3                                                              # number of vectors considered ∈ sizes considered
                cols = findall(!iszero, bvec)
                newstem, type = isf_allsv_from_indices_HnH_r!(Vt, cols, HnH_info, info)         # verification of the subset
            end
    
            # specific technical conditions to update the tree correctly
            if newstem                                                                          # a stem was found
                if type == 0 && p_bvec < p                                                      # for V not Vt, so continued
                    memory = p_bvec                                                             # but some info is stored, this to know where the recursion continued 
                    HnH_info = 1                                                                # now in R^{n+1}
                else    # type == 1                                                             # recursion fully stopped
                    bvec[p_bvec] = 0                                                            # prevent exploring the subtree below
                end
            end
    
            # update the 'binary tree' 
            bvec, p_bvec = isf_binary_tree_update(bvec, p_bvec)
            if bvec[memory] == 0
                memory = p_bvec                                                                 # where the backtracking goes back
                HnH_info = 0                                                                    # recursion comes back to R^n
            end
        end # while p_bvec > 0
    else
        while p_bvec > 0                                                                        # the loop ends when the pointer is zero
            newstem = false
            if 1 < sum(bvec) < r+3                                                              # number of vectors considered ∈ sizes considered
                cols = findall(!iszero, bvec)
                newstem, type = isf_allsv_from_indices_HnH!(Vt, cols, HnH_info, info, options)  # verification of the subset
            end

                                                                                                # specific technical conditions to update the tree correctly
            if newstem                                                                          # a stem was found
                if type == 0 && p_bvec < p                                                      # for V not Vt, so continued
                    memory = p_bvec                                                             # but some info is stored, this to know where the recursion continued 
                    HnH_info = 1                                                                # now in R^{n+1}
                else    # type == 1                                                             # recursion fully stopped
                    bvec[p_bvec] = 0                                                            # prevent exploring the subtree below
                end
            end

            bvec, p_bvec = isf_binary_tree_update(bvec, p_bvec)                                 # update the 'binary tree' 
            if bvec[memory] == 0
                memory = p_bvec                                                                 # where the backtracking goes back
                HnH_info = 0                                                                    # recursion comes back to R^n
            end
        end # while p_bvec > 0
    end

    ## now everything has been checked, purge the lists of their duplicates
    info.stems_sym_init = unique(info.stems_sym_init, dims=2)                                   # purged list
    nb_stems_sym = size(info.stems_sym_init,2)                                                  # actual number of  symmetric stem vectors
    info.nb_duplicated_stems_sym = info.nb_stems_sym - nb_stems_sym
    info.nb_stems_sym = nb_stems_sym

    info.stems_asym_init = unique(info.stems_asym_init, dims=2)                                 # purged list
    nb_stems_asym = size(info.stems_asym_init,2)                                                # actual number of asymmetric stem vectors
    info.nb_duplicated_stems_asym = info.nb_stems_asym - nb_stems_asym
    info.nb_stems_asym = nb_stems_asym

    info.stems_sym = info.stems_sym_init                                                        # the Elastic part will disappear by itself
    info.stems_asym = info.stems_asym_init                                                      # the Elastic part will disappear by itself
    info.stem_sizes_sym  = vec(sum(abs.(info.stems_sym), dims=1))                               # the sum returns a 1 x [size] matrix, not a vector
    info.stem_sizes_asym = vec(sum(abs.(info.stems_asym), dims=1))                              # the sum returns a 1 x [size] matrix, not a vector

    ### If heuristic C (options.bestv=3) is not used, a particular reordering of the stem vectors is used.
    # Said quickly, the idea is that, in this case, the order of the vectors added is fixed by perm.
    # Thus, the indices will be added in the order chosen by perm (thus through the QR factorization).
    # Then, consider a stem vector is given (the permutation forgotten to simplify the notation) by
    # 1-2-3; then one would like to compare a given sign vector s^k = (s_1, ... s_k) with it before
    # comparing with a stem vector that has indices p-2 p-1 p for instance.

    # The first array in the list corresponds to the ones having a 0 in last (taking perm into account)
    # position, the second array the ones with zeros in last two positions and so on.

    if options.bestv == 0 && p > r+1
        info.stem_zero_indices_sym = [Vector{Int64}(undef, 0) for _ in 1:(p-r-1)]
        info.stem_zero_indices_asym = [Vector{Int64}(undef, 0) for _ in 1:(p-r-1)]

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