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

### 
# Multiple (very similar) variants are used for the H, nH, HnH cases.
# The details change a little bit, for technical reasons (see isf_rec_nod_HnH)
# In particular the HnH versions returns a bit more information, 
# used in the main recursive procedure. 
# See also isf_cover_stem.jl for the natural version. 

"""
`b = isf_cover_stem_H_bis(info, options, values)`

Checks whether the sign vector given by svec and perm, in its compacted form, 
covers any of the stem vectors, returning a boolean. 
As the product with the stem matrix is already 'computed', only the 
absolute value and the subtraction are necessary. For more details, 
see `isf_cover_stem.jl` and the functions inside. 

The 'asymmetric' (nH) and 'mixed' (HnH) versions also exist in the same file.
The other version is in `isf_cover_stem.jl`.
"""
function isf_cover_stem_H_bis(svec::Vector{Int64}, info::Info, options::Options)

    p = size(info.stems_sym,1)                  # all are of size p
    ns = length(svec)                           # number of signs in svec

    if options.bestv > 0 || ns == p
        J = Vector(1:info.nb_stems_sym)
    else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        J = info.stem_zero_indices_sym[p-ns]
    end

    # covering is done by the product matrix (of stems) * current vector stored in info.svsprod_sym
    # here with an absolute value the stem vectors are symmetric
    b = 0 in (abs.(info.svsprod_sym[J]) - info.stem_sizes_sym[J])
    return b
end

#-----------------------------------------------------------------------

"""
`b = isf_cover_stem_nH_bis(info, options, values)`

Checks whether the current svec, in its compacted form with the stem
vectors, covers any of them, returning a boolean. 
As the product with the stem matrix is already 'computed', only the 
absolute value and the subtraction are necessary. For more details, 
see `isf_cover_stem.jl` and the functions inside. 
There are two 'current' vector, one for the syms and the other for the asyms. 

The 'symmetric' (H) and 'mixed' (HnH) versions also exist in the same file.
The other version is in `isf_cover_stem.jl`.
"""
function isf_cover_stem_nH_bis(svec::Vector{Int64}, info::Info, options::Options)

    p = size(info.stems_asym,1)                 # all are of size p
    ns = length(svec)                           # number of signs in svec

    if options.bestv > 0 || ns == p
        Ja = Vector(1:info.nb_stems_asym)
    else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        Ja = info.stem_zero_indices_asym[p-ns]
    end

    # covering is done by the product matrix (of stems) * current vector stored in info.svsprod_asym
    # here without absolute value the stem vectors are asymmetric
    b = 0 in (info.svsprod_asym[Ja] - info.stem_sizes_asym[Ja])
    b && return b
    # else => test sym part

    if options.bestv > 0 || ns == p
        Js = Vector(1:info.nb_stems_sym)
    else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        Js = info.stem_zero_indices_sym[p-ns]
    end

    # covering is done by the product matrix (of stems) * current vector stored in info.svsprod_sym
    # here with an absolute value the stem vectors are symmetric
    b = 0 in (abs.(info.svsprod_sym[Js]) - info.stem_sizes_sym[Js])
    return b
end

#-----------------------------------------------------------------------

"""
`b, detail, sign = isf_cover_stem_HnH_bis(HnH_info, info, options, values)`

Checks whether the current svec, in its compacted form with the stem
vectors, covers any of them, returning a boolean and additional information. 
As the product with the stem matrix is already 'computed', only the 
absolute value and the subtraction are necessary. For more details, 
see `isf_cover_stem.jl` and the functions inside. 
There are two 'current' vector, one for the syms and the other for the asyms. 

The 'symmetric' (H) and 'asymmetric' (nH) versions also exist in the same file.
The other version is in `isf_cover_stem.jl`.
"""
function isf_cover_stem_HnH_bis(svec::Vector{Int64}, HnH_info::Int64, info::Info, options::Options)

    p = size(info.stems_sym,1)                      # all are of size p (might be an error if only the ones from somesv... are computed and are all are asymmetric)
    ns = length(svec)

    if options.bestv > 0 || ns == p
        Ja = Vector(1:info.nb_stems_asym)
    else                                            # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        Ja = info.stem_zero_indices_asym[p-ns]
    end

    # covering is done by the product matrix (of stems) * current vector stored in info.svsprod_sym
    # here with an absolute value the stem vectors are used symmetrically
    b = 0 in (abs.(info.svsprod_asym[Ja]) - info.stem_sizes_asym[Ja])
    b && return b, "asym", 0

    ## if a test is successful, it means the current vector covers a stem vector of Vt, 
    ## so the recursion can be (fully) stopped: this is done through the "asym" return
    if HnH_info == 0        # if HnH_info == 1, the recurrence is already in the advanced state so the stems form V (not Vt) are not checked
        if options.bestv > 0 || ns == p
            Js = Vector(1:info.nb_stems_sym)
        else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
            Js = info.stem_zero_indices_sym[p-ns]
        end

        # covering is done by the product matrix (of stems) * current vector stored in info.svsprod_sym
        # here with an absolute value the stem vectors are symmetric
        if !isempty(Js)
            # the index is needed here, so no b = 0 in ...
            j = findfirst(iszero, abs.(info.svsprod_sym[Js]) - info.stem_sizes_sym[Js])
            !isnothing(j) && return true, "sym", sign(info.svsprod_sym[Js[j]])
        end
    end
    return false, "", 0
end