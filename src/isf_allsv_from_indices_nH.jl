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

# 'A priori', there should be a majority of asymmetric 
# stem vectors: they can be 'symmetric' only if they verify 
# an additional equality (thus doesn't happen with random data).

# For the symmetric ones, a convention sign(Z[I[1]])*sign.(Z[I]) is used
# to ensure the first nonzero coordinate is 1: this has no "intrinsic"
# justification but makes the overall list smaller after using unique:
# otherwise there can be stem and -stem in the list, with identical role
# but makes a longer list (not detected by unique).
# This contingency depends on nullspace and cannot be dealt with a priori. 

# For the asymmetric ones, this convention would be an error as stem and
# -stem are not both stem vectors. It depends on the right-hand sides:
# if V*Z = 0, then VS * SZ = 0, but SZ verifies Motzkin's alternative if 
# TS * SZ = TZ >= 0, thus sign(dot(Z, Vt[n+1,cols])) * sign.(Z[I])
# note that this expression is the same if Z <- -Z
# (and the sign is not zero because otherwise it is the symmetric case).

"""
`newstem = isf_allsv_from_indices_nH!(Vt, cols, info, options)`

Modifies info and returns a boolean == true if there is a new stem,
and an additional information for the nature of the stem vector. 
Verifies if the columns 'cols' of V form a circuit (of V(t)'s matroid). 
Computes & checks the dimension of the nullspace of the submatrix without rhs's:

-if different from one (= 0), nothing is done. 
-otherwise, in means these columns contain a circuit: compared to the 
H case, one does not necessarily have +... and -... that are satisfactory. 
Both might exist depending on the values of the right-hand sides, which
is verified afterwards. Then the corresponding stem vector is formed and 
stored, in info.stems_sym or info.stems_asym depending on the right-hand sides. 

Instead of checking if the stem vector exists, it is added and 
duplicates are removed at the end by the 'unique' function. 

The 'symmetric' (H) and 'compact' (HnH) versions are in their respective 
functions (files) `isf_allsv_from_indices_H!(.jl)` and `isf_allsv_from_indices_HnH!(.jl)`.
More technical details can be found at the beginning of the file. 

The computations can be done differently (and often slightly faster), see
`isf_allsv_nH_Wechelon(_rec).jl`.
"""
function isf_allsv_from_indices_nH!(Vt::Matrix, cols::Vector{Int64}, info::Info, options::Options)

    n = size(Vt, 1) - 1                                     # the full matrix is given in argument, we need the 'real' dimension n
    Z = nullspace(Vt[1:n,cols])

    if size(Z,2) != 1                                       # nullity != 1, nothing is done
        return false                                        # the second value is unused in this case
    end

    I = findall(x -> abs(x) > options.tol_coordinates, Z)   # indices of the nonzero components - the real circuit
    stem = zeros(Int, size(Vt,2))
    value = dot(Z, Vt[n+1,cols])                            # quantity determining (a)symmetry
    if abs(value) > options.tol_nonzero_q                   # considered asymmetric so uses the asymmetric convention
        stem[cols[I]] = Int(sign(value))*Int.(sign.(Z[I]))  # the convention is different - one has to take 
        append!(info.stems_asym_init, stem)                 # the appropriate convention, given by this formula
        info.nb_stems_asym += 1                             # updates of info                     
    else                                                    # considered symmetric, uses the symmetric convention
        stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))# which is the convention using the first sign
        append!(info.stems_sym_init, stem)                  # updates of info
        info.nb_stems_sym += 1
    end

    return true
end

# rational version
"""
`newstem = isf_allsv_from_indices_nH_r!(Vt, cols, info)`

Rational version of `isf_allsv_from_indices_nH!`.Uses LinearAlgebraX. 
Modifies info and returns a boolean == true if there is a new stem,
and an additional information for the nature of the stem vector. 
Verifies if the columns 'cols' of V form a circuit (of Vt's matroid). 
Computes & checks the dimension of the nullspace of the submatrix without rhs's:

-if different from one (= 0), nothing is done. 
-otherwise, in means these columns contain a circuit: compared to the 
H case, one does not necessarily have +... and -... that are satisfactory. 
Both might exist depending on the values of the right-hand sides, which
is verified afterwards. Then the corresponding stem vector is formed and 
stored, in info.stems_sym or info.stems_asym depending on the right-hand sides. 

Instead of checking if the stem vector exists, it is added and 
duplicates are removed at the end through the 'unique' function. 

The 'symmetric' (H) and 'compact' (HnH) versions are in their respective 
functions (files) `isf_allsv_from_indices_H_r!(.jl)` and `isf_allsv_from_indices_HnH_r!(.jl)`.
More technical details can be found at the beginning of the file. 

The computations can be done differently (and often slightly faster), see
`isf_allsv_nH_Wechelon(_rec).jl` and its rational version.
"""
function isf_allsv_from_indices_nH_r!(Vt::Matrix, cols::Vector{Int64}, info::Info)

    n = size(Vt, 1) - 1                                     # the full matrix is given in argument, we need the 'real' n dimension
    Z = nullspacex(Vt[1:n,cols]);                           # use of LinearAlgebraX

    if size(Z,2) != 1                                       # nullity != 1, nothing is done
        return false                                        # the second value is unused in this case
    end

    I = findall(x -> x != 0, Z)                             # indices of the nonzero components - the real circuit
    stem = zeros(Int, size(Vt,2))
    value = dot(Z, Vt[n+1,cols])                            # quantity determining (a)symmetry
    if abs(value) > 0                                       # considered asymmetric so uses the asymmetric convention
        stem[cols[I]] = Int(sign(value))*Int.(sign.(Z[I]))  # the convention is different - one has to take 
        append!(info.stems_asym_init, stem)                 # the appropriate convention, given by this formula
        info.nb_stems_asym += 1                             # updates of info                   
    else                                                    # considered symmetric, uses the symmetric convention
        stem[cols[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))# which is the convention using the first sign
        append!(info.stems_sym_init, stem)                  # updates of info
        info.nb_stems_sym += 1
    end

    return true
end