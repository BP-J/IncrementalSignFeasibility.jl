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
`isf_allsv_H_Wechelon_rec!(V, perm, cols, starting_echelon, info, options)`

Inspired from the computations of TOPCOM with the echelonned form. 

Recursive call on the subset size to possibly modify info. 
In each 'first' call, the size is 2 so it checks a subset of size 3. 
The echelon form is used purely on the submatrices, and nullspace is called
when the last line becomes zero. 

The 'asymmetric' (nH) and 'compact' (HnH) versions are in their respective 
functions (files) `isf_allsv_nH_Wechelon_rec!(.jl)`, `isf_allsv_HnH_Wechelon_rec!(.jl)`.

For additional details, see `isf_allsv_from_indices_H(.jl)`, which is 
the more natural algorithm, but often slightly slower. 
"""
function isf_allsv_H_Wechelon_rec!(V::Matrix, cols::Vector{Int64}, perm::Vector{Int64}, starting_echelon::Matrix, info::Info, options::Options)

    n, p = size(V)  
    for new_index in maximum(cols)+1:p                                              # launch the following recursive calls
        cols_plus = [cols;new_index]
        temp_echelon_matrix = echelonned_form([starting_echelon ; V[:,new_index]']) # echelonned matrix with the new vector

        if norm(temp_echelon_matrix[length(cols)+1,1:n]) < options.tol_nonzero_q    # considered to be zero, so there is a circuit/stem to find in the current subset
            Z = nullspace(V[:,cols_plus]);
            I = findall(x -> abs(x) > options.tol_coordinates, Z)                   # indices of the nonzero components - the real circuit
            stem = zeros(Int, p)
            stem[cols_plus[I]] = sign(Z[I[1]])*sign.(Z[I])                          # construction of the stem vector with convention
            append!(info.stems_sym_init, stem)                                      # update of info
            info.nb_stems_sym += 1
        else                                                                        # continue the recursion, subset still independent even with the new one
            isf_allsv_H_Wechelon_rec!(V, cols_plus, perm, temp_echelon_matrix, info, options)
        end
    end
end

# rational version
"""
`isf_allsv_H_Wechelon_rec!(V, perm, cols, starting_echelon, info, options)`

Inspired from the computations of TOPCOM with the echelonned form. 

Rational version of `isf_allsv_from_indices_H!`. Uses LinearAlgebraX. 
Recursive call on the subset size to eventually modify info. 
In each 'first' call, the size is 2 so it checks a subset of size 3. 
The echelon form is used purely on the submatrices, and nullspace is called
when the last line becomes zero. 

The 'asymmetric' (nH) and 'compact' (HnH) versions are in their respective 
functions (files) `isf_allsv_nH_Wechelon_rec!(.jl)`, `isf_allsv_HnH_Wechelon_rec!(.jl)`.

For additional details, see `isf_allsv_from_indices_H(.jl)`, which is 
the more natural algorithm, but often slightly slower. 
"""
function isf_allsv_H_Wechelon_rec_r!(V::Matrix, cols::Vector{Int64}, perm::Vector{Int64}, starting_echelon::Matrix, info::Info, options::Options)

    n, p = size(V)
    s = length(cols)        
    for new_index in maximum(cols)+1:p                                              # launch the following recursive calls
        cols_plus = [cols;new_index]
        temp_echelon_matrix = echelonned_form([starting_echelon ; V[:,new_index]']) # echelonned matrix with the new vector

        if norm(temp_echelon_matrix[s+1,1:n]) == 0                                  # considered to be zero, so there is a circuit/stem to find in the current subset
            Z = nullspacex(V[:,cols_plus]);                                         # use of LinearAlgebraX
            I = findall(x -> x != 0, Z)                                             # indices of the nonzero components - the real circuit
            stem = zeros(Int, p)
            stem[cols_plus[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))               # construction of the stem vector with convention
            append!(info.stems_sym_init, stem)                                      # update of info
            info.nb_stems_sym += 1
        else                                                                        # continue the recursion, subset still independent even with the new one
            isf_allsv_H_Wechelon_rec_r!(V, cols_plus, perm, temp_echelon_matrix, info, options)
        end
    end
end