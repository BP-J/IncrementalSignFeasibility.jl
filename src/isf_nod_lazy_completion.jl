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
`isf_nod_completion_H!(stem, base, size, p, info)`

Computes a part of the complementary set from a stem vector.
The function manipulates Booleans (0-1), and distinguishes 
if the first component is in the stem vector: if yes, 
then one fills the gaps in the 2^{p-l} vectors. 
If no, since by symmetry one wants the first component to be 1, 
then one must use stem on half the completed vectors and -stem
on the other half. 

Duplicates the stem vector and its opposite through 'repeat'
then completes the empty parts. Used for the symmetric setting, 
generating only a part of the complementary set.
"""
function isf_nod_completion_H!(stem::SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, base::Vector{Int64}, p::Int64, info::Info)

    stembase = Bool.((stem[base] .+ 1)/2)               # the common thing to all sign vectors covering the stem vector
    M_base = Bool.([1 0])
    comp_ind = setdiff(1:p,base)                        # the other indices
    size = length(base)
    M = ones(Bool, p, 2^(p-size))                       # the matrix of complements
    if 1 in base                                        # completing with 0s and 1s at the indices of comp_ind
        for i in 1:p-size
            M[comp_ind[i],:] = repeat(M_base, inner=(1,2^(p-size-i)), outer=(1,2^(i-1)))
        end
        M[base,:] = repeat(stembase, outer=(1,2^(p-size)))

        union!(info.sclst, eachcol(M))                  # operation on a Set
    else                                                # one must use the stem and its opposite to guarantee the first component remains 1 (for symmetry)
        for i in 1:p-size-1
            M[comp_ind[i+1],:] = repeat(M_base, inner=(1,2^(p-size-i-1)), outer=(1, 2^i))
        end
        M[base,:] = repeat([stembase .~stembase], inner=(1,2^(p-size-1)))

        union!(info.sclst, eachcol(M))                  # operation on a Set 
    end
end

"""
`isf_nod_completion_nH_sym!(stem, base, size, p, info)`

Computes a part of the complementary set from a stem vector. 
Duplicates the stem vector through 'repeat' then completes 
the empty parts. Used for the asymmetric setting, generating 
every sign vector covered by the stem vector and its opposite.
"""
function isf_nod_completion_nH_sym!(stem::SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, base::Vector{Int64}, p::Int64, info::Info)

    stembase = Bool.((stem[base] .+ 1)/2)               # the common thing to all sign vectors covering the stem vector
    M_base = Bool.([1 0])
    comp_ind = setdiff(1:p, base)                       # the other indices
    size = length(base)
    M = ones(Int64, p, 2^(p-size+1))                    # the matrix of complements, of size p x 2^(p-size+1) since its a symmetric stem 
    for i in 1:p-size                                   # completing with 0s and 1s at the indices of comp_ind, using also the stem's opposite 
        M[comp_ind[i],:] = repeat(M_base, inner=(1,2^(p-size-i)), outer=(1,2^i))
    end
    M[base,:] = repeat([stembase .~stembase], inner=(1, 2^(p-size)))

    union!(info.sclst, eachcol(M))                      # operation on a Set 
end


"""
`isf_nod_completion_nH_asym!(stem, base, size, p, info)`

Computes a part of the complementary set from a stem vector.
Duplicates the stem vector through 'repeat'
then completes the empty parts. Used for the asymmetric setting, 
generating every sign vector covered by the stem vector.
"""
function isf_nod_completion_nH_asym!(stem::SubArray{Int64, 1, Matrix{Int64}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}, base::Vector{Int64}, p::Int64, info::Info)

    stembase = Bool.((stem[base] .+ 1)/2)               # the common thing to all sign vectors covering the stem vector
    M_base = Bool.([1 0])
    comp_ind = setdiff(1:p, base)                       # the other indices
    size = length(base)
    M = ones(Int64, p, 2^(p-size))                      # the matrix of complements of size p x 2^(p-size)
    for i in 1:p-size                                   # completing with 0s and 1s at the indices of comp_ind, using also the stem's opposite 
        M[comp_ind[i],:] = repeat(M_base, inner=(1,2^(p-size-i)), outer=(1,2^(i-1)))
    end
    M[base,:] = repeat(stembase, inner=(1, 2^(p-size)))

    union!(info.sclst, eachcol(stemmed))                # operation on a Set 
end

# There is no point in having a HnH variant since no tree is used. 