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
# See also isf_cover_stem_bis.jl for a sometimes faster version. 

"""
`b = isf_cover_stem_H(svec, perm, info, options, values)`

Checks whether the sign vector given by svec and perm covers an element
in info.stems_sym, returning a boolean. 

The stem vectors are in {-1, 0, +1}^p.
The ± 1 mean the indices define the stem and the indices with 0 are unused.
A sign vector s covers a stem vector stem if stem[NonZero] = ± s[NonZero], with 
NonZero being the indices such that stem(i) != 0.

The 'asymmetric' (nH) and 'mixed' (HnH) versions also exist in the same file.
The other version is in `isf_cover_stem_bis.jl`.
"""
function isf_cover_stem_H(svec::Vector{Int64}, perm::Vector{Int64}, info::Info, options::Options)
    
    p = size(info.stems_sym,1)                  # all are of size p (possibly with some 0 components)
    ns = length(svec)                           # number of signs in svec

    if options.bestv > 0 || ns == p
        J = Vector(1:info.nb_stems_sym)
    else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        J = info.stem_zero_indices_sym[p-ns]
    end

    s = zeros(Int,p)
    s[perm[1:ns]] = svec

    # covering is done by the product matrix (of stems) * current vector
    # here with an absolute value the stem vectors are symmetric
    b = 0 in (abs.(info.stems_sym[:,J]' * s) - info.stem_sizes_sym[J])
    return b
end

#-----------------------------------------------------------------------

"""
`b = isf_cover_stem_nH(svec, perm, info, options, values)`

Checks whether the sign vector given by svec and perm covers an element
in info.stems_sym or in info.stems_asym, returning a boolean. 

The stem vectors are in {-1, 0, +1}^p.
The ± 1 mean the indices define the stem and the indices with 0 are unused.
A sign vector s covers a   symmetric stem vector stem if stem[NonZero] = ±s[NonZero], 
with NonZero being the indices such that stem(i) != 0. 
A sign vector s covers an asymmetric stem vector stem if stem[NonZero] = s[NonZero], 
with NonZero being the indices such that stem(i) != 0. 

The 'symmetric' (H) and 'mixed' (HnH) versions also exist in the same file.
The other version is in `isf_cover_stem_bis.jl`.
"""
function isf_cover_stem_nH(svec::Vector{Int64}, perm::Vector{Int64}, info::Info, options::Options)

    p = size(info.stems_asym,1)                 # all are of size p
    ns = length(svec)                           # number of signs in svec

    if options.bestv > 0 || ns == p
        Ja = Vector(1:info.nb_stems_asym)
    else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        Ja = info.stem_zero_indices_asym[p-ns]
    end

    s = zeros(Int, p)
    s[perm[1:ns]] = svec

    # covering is done by the product matrix (of stems) * current vector
    # here without absolute value as the stem vectors are asymmetric
    b = 0 in (info.stems_asym[:,Ja]' * s - info.stem_sizes_asym[Ja])
    b && return b
    # else => test sym part

    if options.bestv > 0 || ns == p
        Js = Vector(1:info.nb_stems_sym)
    else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        Js = info.stem_zero_indices_sym[p-ns]
    end

    # covering is done by the product matrix (of stems) * current vector
    # here with an absolute value the stem vectors are symmetric
    b = 0 in (abs.(info.stems_sym[:,Js]' * s) - info.stem_sizes_sym[Js])
    
    return b
end

#-----------------------------------------------------------------------

"""
`b, detail, sign = isf_cover_stem_HnH(svec, perm, HnH_info, info, options, values)`

Checks whether the sign vector given by svec and perm covers an element
in info.stems_sym or in info.stems_asym, returning a boolean and additional information. 

The stem vectors are in {-1, 0, +1}^p.
The ± 1 mean the indices define the stem and the indices with 0 are unused.
Both types of stem vectors are used "symmetrically" despite their name: 
s covers a   symmetric stem vector stem if stem[NonZero] = ± s[NonZero], with 
NonZero being the indices such that stem(i) != 0. 
s covers an asymmetric stem vector stem if stem[NonZero] = ± s[NonZero], with
NonZero being the indices such that stem(i) != 0. 

The 'symmetric' (H) and 'asymmetric' (nH) versions also exist in the same file.
The other version is in `isf_cover_stem_bis.jl`.
"""
function isf_cover_stem_HnH(svec::Vector{Int64}, perm::Vector{Int64}, HnH_info::Int64, info::Info, options::Options)

    p = size(info.stems_sym,1)                      # all are of size p (might be an error if only the ones from somesv... are computed and are all are asymmetric)
    ns = length(svec)

    s = zeros(Int, p)
    s[perm[1:ns]] = svec

    if options.bestv > 0 || ns == p
        Ja = Vector(1:info.nb_stems_asym)
    else                                            # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
        Ja = info.stem_zero_indices_asym[p-ns]
    end

    # covering is done by the product matrix (of stems) * current vector
    # here with an absolute value the stem vectors are used symmetrically
    b = 0 in abs.(info.stems_asym[:,Ja]' * s) - info.stem_sizes_asym[Ja]
    b && return b, "asym", 0
    
    ## if a test is successful, it means the current vector covers a stem vector of Vt, 
    ## so the recursion can be (fully) stopped: this is done through the "asym" return
    if HnH_info == 0                                # if HnH_info == 1, the recurrence is already in the advanced state so the stems form V (not Vt) are not checked
        if options.bestv > 0 || ns == p
            Js = Vector(1:info.nb_stems_sym)
        else                                        # options.bestv = 0 so we use only some stem vectors (who have appropriate coordinates)
            Js = info.stem_zero_indices_sym[p-ns]
        end

        # covering is done by the product matrix (of stems) * current vector
        # here with an absolute value the stem vectors are used symmetrically
        if !(isempty(Js))
            # the index is needed here, so no b = 0 in ...
            j = findfirst(iszero, abs.(info.stems_sym[:,Js]' * s) - info.stem_sizes_sym[Js])
            !(isnothing(j)) && return true, "sym", sign(s' * info.stems_sym[:,Js[j]])
        end
    end
    return false, "", 0
end