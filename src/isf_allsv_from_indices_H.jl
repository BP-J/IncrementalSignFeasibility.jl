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

### Details
# Instead of checking duplicates at every addition, the code checks 
# at the end for duplicates. For instance, if one is found by
# 3-4-5, then subsets 1-2-3-4-5 and 2-3-4-5 will also find it. 

# info.stems_asym stores the asymmetric stem vectors and 
# is used only in the nH and HnH settings.

# For the symmetric ones, a convention sign(Z[I[1]])*sign.(Z[I]) is used
# to ensure the first nonzero coordinate is 1: this has no "intrinsic"
# justification but makes the overall list smaller after using unique:
# otherwise there can be stem and -stem in the list, with identical role
# but makes a longer list (not detected by unique).
# This contingency depends on nullspace and cannot be dealt with a priori. 

"""
`newstem = isf_allsv_from_indices_H!(V, cols, info, options)`

Modifies info and returns a boolean == true if there is a new stem
that is a subset of 'cols' (a circuit of V's matroid). 
Computes & checks the dimension of the nullspace of the submatrix:

-if different from one (= 0), nothing is done. 
-otherwise, it means these columns contain a circuit. Then, 
the corresponding stem vector is extracted, formed and stored in
info.stems_sym, with the 'first nonzero term = 1' convention, to
avoid having a stem and its opposite (which have the same role).

Instead of checking if the stem vector exists, it is added and 
duplicates are removed at the end by the 'unique' function. 

The 'asymmetric' (nH) and 'compact' (HnH) versions are in their respective
functions (files) `isf_allsv_from_indices_nH!(.jl)` and `isf_allsv_from_indices_HnH(.jl)`. 
More technical details can be found at the beginning of the file. 

The computations can be done differently (and often slightly faster), see
`isf_allsv_H_Wechelon(_rec).jl`.
"""
function isf_allsv_from_indices_H!(V::Matrix, cols::Vector{Int64}, info::Info, options::Options)
    
    Z = nullspace(V[:,cols])
    
    if size(Z,2) != 1                                       # nullity != 1, nothing is done
        return false
    end

    I = findall(x -> abs(x) > options.tol_coordinates, Z)   # indices of the nonzero components - the real circuit
    stem = zeros(Int, size(V,2))
    stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))    # construction of the stem vector with convention
    append!(info.stems_sym_init, stem)                      # updates of info
    info.nb_stems_sym += 1
    return true
end

# rational version
"""
`newstem = isf_allsv_from_indices_H_r!(V, cols, info)`

Rational version of `isf_allsv_from_indices_H!`. Uses LinearAlgebraX. 
Modifies info and returns a boolean == true if there is a new stem
that is a subset of 'cols' (a circuit of V's matroid). 
Computes & checks the dimension of the nullspace of the submatrix:

-if different from one (= 0), nothing is done. 
-otherwise, it means these columns contain a circuit. Then, 
the corresponding stem vector is extracted, formed and stored in
info.stems_sym, with the 'first nonzero term = 1' convention, to
avoid having a stem and its opposite (which have the same role).

Instead of checking if the stem vector exists, it is added and 
duplicates are removed at the end through the 'unique' function. 

The 'asymmetric' (nH) and 'compact' (HnH) versions are in their respective
functions (files) `isf_allsv_from_indices_nH_r!(.jl)` and `isf_allsv_from_indices_HnH_r!(.jl)`. 
More technical details can be found at the beginning of the file. 

The computations can be done differently (and often slightly faster), see
`isf_allsv_H_Wechelon(_rec).jl` and its rational version.
"""
function isf_allsv_from_indices_H_r!(V::Matrix, cols::Vector{Int64}, info::Info)
    
    Z = nullspacex(V[:,cols]);                              # use of LinearAlgebraX
    
    if size(Z,2) != 1                                       # nullity != 1, nothing is done
        return false
    end

    I = findall(x -> x != 0, Z)                             # indices of the nonzero components - the real circuit
    stem = zeros(Int, size(V,2))
    stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))    # construction of the stem vector with convention
    append!(info.stems_sym_init, stem)                      # updates of info
    info.nb_stems_sym += 1
    return true
end