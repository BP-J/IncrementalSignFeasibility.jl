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

# There are many symmetric and asymmetric stem vectors in this case.
# !!! The convention is different with the paper cited above !!!

# For the "symmetric" stem vectors that are the stem vectors of V
# (and not of Vt), another convention is used: the formula 
# [stem] = sign(dot(Z, Vt[n+1,cols])) * sign.(Z[I]) where Z=nullspace(V)
# (V[:,cols] noted V for simplicity) comes from the following: 
# if the columns cols have nullity one, then the null space is given by
# Vect(Z), and the stem vector is given by sign.(Z) 
# (for coordinates 'large enough' in absolute value); 
# then V*Z = 0, but for a later use it is better to have the "real" sign vector
# i.e., the one given in the non-homogeneous case, i.e., oriented by 
# the sign of the scalar product in the formula. This convention, 
# different from the other one, also prevents from having stem and -stem.
# (Basically Z must appear twice in the formula so that Z and -Z returns the same result).

# For the "asymmetric" ones, i.e., the stem vectors of Vt seen in R^{n+1}, 
# they are used in a symmetric way so the orientation does not matter, 
# and the usual sign(Z[I[1]])*sign.(Z[I]) 'first non-zero element is +1' 
# convention is used (the symmetric convention of the H and nH cases).

# This computation can be done differently, for instance by using
# 2 null spaces computations, though the current seemed faster.

"""
`newstem, type = isf_allsv_from_indices_HnH!(Vt, cols, info, options)`

Modifies info and returns a boolean == true if there is a new stem; 
type is an additional information for the nature of the stem vector. 
Verifies if the columns 'cols' of Vt form a circuit (of Vt's matroid). 
Computes & checks the dimension of the nullspace of the submatrix.

First, if a subset of the columns is a circuit in R^n already, 
then it is searched in R^{n+1}.

If searched in R^n, when the null space is of size 1, the corresponding stem
vector is formed and stored, possibly already in R^{n+1} if the rhs's coincide. 
Otherwise a special sign convention is used to determine the nature of the 
current subbranch. If searched in R^{n+1}, when the null space is of size 1, it is 
necessarily an 'asymmetric' stem vector as the right-hand sides are already considered. 

Instead of checking if the stem vector exists, it is added and 
duplicates are removed at the end by the 'unique' function. 

The 'symmetric' (H) and 'asymmetric' (nH) versions are in their respective 
functions (files) `isf_allsv_from_indices_H!(.jl)` and `isf_allsv_from_indices_nH!(.jl)`.
More technical details can be found at the beginning of the file. 

The computations can be done differently (and often slightly faster), see 
`isf_allsv_HnH_Wechelon(_rec).jl`. 
"""
function isf_allsv_from_indices_HnH!(Vt::Matrix, cols::Vector{Int64}, HnH_info::Int64, info::Info, options::Options)

    p = size(Vt,2)
    n = size(Vt, 1) - 1 + HnH_info                                  # the dimension depends on the recurrence state

    Z = nullspace(Vt[1:n,cols])                                     # the rhs as new coordinates are checked later
    S = size(Z,2)

    if HnH_info == 0                                                # computation still in R^{n}
        if S == 0                                                   # dimension too small, moving to the next subset
            return false, 0
        elseif S == 1
            newstem = true
            I = findall(x -> abs(x) > options.tol_coordinates, Z)   # indices of the nonzero components - the real circuit
            stem = zeros(Int, p)
            value = dot(Vt[n+1,cols], Z)                            # quantity determining (a)symmetry
            if abs(value) <= options.tol_nonzero_q                  # considered in Vt's nullspace
                stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))# one could also use the sign(T*...) convention but not necessary
                append!(info.stems_asym_init, stem)                 # updates of info
                info.nb_stems_asym += 1
                type = 1                                            # "asymmetric" stem vector
            else
                stem[cols[I]] = Int(sign(value))*Int.(sign.(Z[I]))  # sign convention, see above
                # the scalar product cannot be close to 0, otherwise Z would be in the 
                # nullspace of Vt which is the case just above
                append!(info.stems_sym_init, stem)                  # updates of info
                info.nb_stems_sym += 1
                type = 0                                            # "symmetric" stem vector
            end
        else                                                        # probably unnecessary
            return false, 0
        end
    else                                                            # computation in R^{n+1}, already detected something in R^n before
        if S == 0                                                   # dimension too small, moving to the next one
            return false, 1                                         # here type is still 1
        elseif S == 1
            newstem = true
            I = findall(x -> abs(x) > options.tol_coordinates, Z)
            stem = zeros(Int, p)
            stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))    # sign convention, see above
            append!(info.stems_asym_init, stem)                     # updates of info
            info.nb_stems_asym += 1
            type = 1
        else                                                        # probably unnecessary
            return false, 1
        end
    end
    return newstem, type
end

# rational version
"""
`newstem, type = isf_allsv_from_indices_HnH_r!(Vt, cols, info)`

Rational version of `isf_allsv_from_indices_HnH!`. Uses LinearAlgebraX. 
Modifies info and returns a boolean == true if there is a new stem; 
type is an additional information for the nature of the stem vector. 
Verifies if the columns 'cols' of Vt form a circuit (of Vt's matroid). 
Computes & checks the dimension of the nullspace of the submatrix.

First, if a subset of the columns is a circuit in R^n already, 
then it is searched in R^{n+1}.

If searched in R^n, when the null space is of size 1, the corresponding stem
vector is formed and stored, possibly already in R^{n+1} if the rhs's coincide. 
Otherwise a special sign convention is used to determine the nature of the 
current subbranch. If searched in R^{n+1}, when the null space is of size 1, it is 
necessarily an 'asymmetric' stem vector as the right-hand sides are already considered. 

Instead of checking if the stem vector exists, it is added and 
duplicates are removed at the end through the 'unique' function. 

The 'symmetric' (H) and 'asymmetric' (nH) versions are in their respective 
functions (files) `isf_allsv_from_indices_H_r!(.jl)` and `isf_allsv_from_indices_nH_r!(.jl)`.
More technical details can be found at the beginning of the file. 

The computations can be done differently (and often slightly faster), see 
`isf_allsv_HnH_Wechelon(_rec).jl` and its rational version.
"""
function isf_allsv_from_indices_HnH_r!(Vt::Matrix, cols::Vector{Int64}, HnH_info::Int64, info::Info)

    p = size(Vt,2)
    n = size(Vt, 1) - 1 + HnH_info                                  # the dimension depends on the recurrence state

    Z = nullspacex(Vt[1:n,cols])                                    # use of LinearAlgebraX; the rhs as new coordinates are checked later
    S = size(Z,2)

    if HnH_info == 0                                                # computation still in R^{n}
        if S == 0                                                   # dimension too small, moving to the next subset
            return false, 0
        elseif S == 1
            newstem = true
            I = findall(x -> x != 0, Z)                             # indices of the nonzero components - the real circuit
            stem = zeros(Int, p)
            value = dot(Vt[n+1,cols], Z)                            # quantity determining (a)symmetry
            if abs(value) == 0                                      # considered in Vt's nullspace
                stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))# one could also use the sign(T*...) convention but not necessary
                append!(info.stems_asym_init, stem)                 # updates of info
                info.nb_stems_asym += 1
                type = 1                                            # "asymmetric" stem vector
            else
                stem[cols[I]] = Int(sign(value))*Int.(sign.(Z[I]))  # sign convention, see above
                # the scalar product cannot be close to 0, otherwise Z would be in the 
                # nullspace of Vt which is the case just above
                append!(info.stems_sym_init, stem)                  # updates of info
                info.nb_stems_sym += 1
                type = 0                                            # "symmetric" stem vector
            end
        else                                                        # probably unnecessary
            return false, 0
        end
    else                                                            # computation in R^{n+1}, already detected something in R^n before
        if S == 0                                                   # dimension too small, moving to the next one
            return false, 1                                         # here type is still 1
        elseif S == 1
            newstem = true
            I = findall(x -> x != 0, Z)
            stem = zeros(Int, p)
            stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))    # sign convention, see above
            append!(info.stems_asym_init, stem)                     # updates of info
            info.nb_stems_asym += 1
            type = 1
        else                                                        # probably unnecessary
            return false, 1
        end
    end
    return newstem, type
end