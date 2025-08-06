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
`isf_allsv_HnH_Wechelon!(Vt, perm, info, options, values)`

Prepares the recursive calls to `isf_allsv_HnH_Wechelon_rec!`
Modifies info through a recursion on the size of the current subset.

More precisions are found in `isf_allsv_HnH_Wechelon_rec.jl`.
The 'symmetric' (H) and 'asymmetric' (nH) versions are in their respective 
functions (files) `isf_allsv_H_Wechelon!(.jl)`, `isf_allsv_nH_Wechelon!(.jl)`.

For additional details, see `isf_allsv_HnH(.jl)`, which is a more natural 
introduction, but often slightly slower.
"""
function isf_allsv_HnH_Wechelon!(Vt::Matrix, perm::Vector{Int64}, info::Info, options::Options, values::Values)

    p = size(Vt, 2)
    n = size(Vt, 1) - 1
    if p < 3
        println("isf_allsv_HnH_Wechelon!: function is not made for small values of p = $(p)")
        info.flag = values.fail_technicality
    end
    r = info.r

# ----------------------------------------------------------------------------------------------------
# Compute all the stem vector of a matrix without permuting the columns (version with a binary tree)
# ----------------------------------------------------------------------------------------------------

    # computes the binomial(p, 1) starts, then calls the recursive process on each one
    # the affine case can consider pairwise stems, as the right-hand sides intervene
    if options.rational                                                 # test before to avoid doing it every time
        for i in 1:p
            isf_allsv_HnH_Wechelon_rec_r!(Vt, [i], perm, echelonned_form(Vt[1:n,i]'), true, info, options)
        end
    else
        for i in 1:p
            isf_allsv_HnH_Wechelon_rec!(Vt, [i], perm, echelonned_form(Vt[1:n,i]'), true, info, options)
        end
    end
                
    ## now everything has been checked, purge the lists of their duplicates
    info.stems_sym_init = unique(info.stems_sym_init, dims=2)           # purged list
    nb_stems_sym = size(info.stems_sym_init,2)                          # actual number of  symmetric stem vectors
    info.nb_duplicated_stems_sym = info.nb_stems_sym - nb_stems_sym
    info.nb_stems_sym = nb_stems_sym

    info.stems_asym_init = unique(info.stems_asym_init, dims=2)         # purged list
    nb_stems_asym = size(info.stems_asym_init,2)                        # actual number of asymmetric stem vectors
    info.nb_duplicated_stems_asym = info.nb_stems_asym - nb_stems_asym
    info.nb_stems_asym = nb_stems_asym

    info.stems_sym = info.stems_sym_init                                # the Elastic part will disappear by itself
    info.stems_asym = info.stems_asym_init                              # the Elastic part will disappear by itself
    info.stem_sizes_sym  = vec(sum(abs.(info.stems_sym), dims=1))       # the sum returns a 1 x [size] matrix, not a vector
    info.stem_sizes_asym = vec(sum(abs.(info.stems_asym), dims=1))      # the sum returns a 1 x [size] matrix, not a vector

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