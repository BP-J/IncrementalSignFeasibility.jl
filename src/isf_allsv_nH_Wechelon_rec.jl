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
`isf_allsv_nH_Wechelon_rec!(Vt, perm, cols, starting_echelon, info, options)`

Inspired from the computations of TOPCOM with the echelonned form.

Recursive call on the subset size to possibly modify info. 
In each 'first' call, the size is 1 so it checks a subset of size 2. 
The echelon form is used purely on the submatrices, and nullspace is called
when the last line becomes zero. 
There is a verification whether the stem vector might be symmetric (unlikely), 
which uses the right-hand sides. 

The 'symmetric' (H) and 'compact' (HnH) versions are in their respective 
functions (files) `isf_allsv_H_Wechelon_rec!(.jl)`, `isf_allsv_HnH_Wechelon_rec!(.jl)`.

For additional details, see `isf_allsv_from_indices_nH(.jl)`, which is 
the more natural algorithm, but often slightly slower. 
"""
function isf_allsv_nH_Wechelon_rec!(Vt::Matrix, cols::Vector{Int64}, perm::Vector{Int64}, starting_echelon::Matrix, info::Info, options::Options)

    p = size(Vt, 2)
    n = size(Vt, 1) - 1
    for new_index in maximum(cols)+1:p                                                  # launch the following recursive calls
        cols_plus = [cols;new_index]
        temp_echelon_matrix = echelonned_form([starting_echelon ; Vt[1:n,new_index]' ]) # echelonned matrix with the new vector
        
        if norm(temp_echelon_matrix[length(cols)+1,1:n]) < options.tol_nonzero_q        # considered to be zero, so there is a circuit/stem to find in the current subset
            Z = nullspace(Vt[1:n,cols_plus])
            I = findall(x -> abs(x) > options.tol_coordinates, Z)                       # indices of the nonzero components - the real circuit
            stem = zeros(Int, size(Vt,2))
            value = dot(Z, Vt[n+1,cols_plus])                                           # quantity determining (a)symmetry
            if abs(value) > options.tol_nonzero_q                                       # considered asymmetric so uses the asymmetric convention
                stem[cols_plus[I]] = sign(value)*sign.(Z[I])                            # the convention is different - one has to take 
                append!(info.stems_asym_init, stem)                                     # the right convention, given by this formula
                info.nb_stems_asym += 1                                                 # update of info
            else                                                                        # considered that both stem vectors exist 
                stem[cols_plus[I]] = sign(Z[I[1]])*sign.(Z[I])                          # construction of the stem vector with convention
                append!(info.stems_sym_init, stem)                                      # update of info
                info.nb_stems_sym += 1
            end
        else                                                                            # continue the recursion, subset still independent even with the new one
            isf_allsv_nH_Wechelon_rec!(Vt, cols_plus, perm, temp_echelon_matrix, info, options)
        end
    end
end

# rational version
"""
`isf_allsv_nH_Wechelon_rec!(Vt, perm, cols, starting_echelon, info, options)`

Inspired from the computations of TOPCOM with the echelonned form.

Rational version of `isf_allsv_from_indices_nH!`. Uses LinearAlgebraX. 
Recursive call on the subset size to possibly modify info. 
In each 'first' call, the size is 1 so it checks a subset of size 2. 
The echelon form is used purely on the submatrices, and nullspace is called
when the last line becomes zero. 
There is a verification whether the stem vector might be symmetric (unlikely), 
which uses the right-hand sides. 

The 'symmetric' (H) and 'compact' (HnH) versions are in their respective 
functions (files) `isf_allsv_H_Wechelon_rec!(.jl)`, `isf_allsv_HnH_Wechelon_rec!(.jl)`.

For additional details, see `isf_allsv_from_indices_nH(.jl)`, which is 
the more natural algorithm, but often slightly slower. 
"""
function isf_allsv_nH_Wechelon_rec_r!(Vt::Matrix, cols::Vector{Int64}, perm::Vector{Int64}, starting_echelon::Matrix, info::Info, options::Options)

    p = size(Vt, 2)
    n = size(Vt, 1) - 1
    s = length(cols)
    for new_index in maximum(cols)+1:p                                                  # launch the following recursive calls
        cols_plus = [cols;new_index]
        temp_echelon_matrix = echelonned_form([starting_echelon ; Vt[1:n,new_index]'])  # echelonned matrix with the new vector

        if norm(temp_echelon_matrix[s+1,1:n]) == 0                                      # so there is a circuit/stem to find in the current subset
            Z = nullspacex(Vt[1:n,cols_plus])                                           # use of LinearAlgebraX
            I = findall(x -> x != 0, Z)                                                 # indices of the nonzero components - the real circuit
            stem = zeros(Int, size(Vt,2))
            value = dot(Z, Vt[n+1,cols_plus])                                           # quantity determining (a)symmetry
            if abs(value) > 0                                                           # considered asymmetric so uses the asymmetric convention
                stem[cols_plus[I]] = Int(sign(value))*Int.(sign.(Z[I]))                 # the convention is different - one has to take 
                append!(info.stems_asym_init, stem)                                     # the right convention, given by this formula
                info.nb_stems_asym += 1                                                 # update of info
            else                                                                        # considered that both stem vectors exist 
                stem[cols_plus[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))               # construction of the stem vector with convention
                append!(info.stems_sym_init, stem)                                      # update of info
                info.nb_stems_sym += 1
            end
        else                                                                            # continue the recursion, subset still independent even with the new one
            isf_allsv_nH_Wechelon_rec_r!(Vt, cols_plus, perm, temp_echelon_matrix, info, options)
        end
    end
end