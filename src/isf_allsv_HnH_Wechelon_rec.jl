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
`isf_allsv_HnH_Wechelon_rec!(Vt, perm, cols, starting_echelon, info, options)`

Inspired from the computations of TOPCOM with the echelonned form. 

Recursive call on the subset size to possibly modify info. 
In each 'first' call, the size is 1 so it checks a subset of size 2. 
The echelon form is used purely on the submatrices, and nullspace is called
when the last line becomes zero. 
While the current recursion is in R^n, it is similar as in the nH case. 
When it is in R^{n+1}, the last line is already zero and the 
stem vector found wasn't one in R^{n+1} (otherwise it is already known). 
Thus when the _two_ last lines are zero, in R^n because this is where 
the echelon form is computed, necessarily in R^{n+1} there is a circuit,
which is obtained through the combination of both. 

The 'symmetric' (H) and 'asymmetric' (nH) versions are in their respective 
functions (files) `isf_allsv_H_Wechelon_rec!(.jl)`, `isf_allsv_nH_Wechelon_rec!(.jl)`.

For additional details, see `isf_allsv_from_indices_HnH(.jl)`, which is a more natural 
introduction, but often slightly slower. 
"""
function isf_allsv_HnH_Wechelon_rec!(Vt::Matrix, cols::Vector{Int64}, perm::Vector{Int64}, starting_echelon::Matrix, memory::Bool, info::Info, options::Options)

    p = size(Vt, 2)
    n = size(Vt, 1) - 1

    for new_index in maximum(cols)+1:p                                                  # launch the following recursive calls
        cols_plus = [cols;new_index]
        temp_echelon_matrix = echelonned_form([starting_echelon ; Vt[1:n,new_index]'])  # echelonned matrix with the new vector

        if memory                                                                       # means previous vectors were still independent (natural phase, R^n)
            if norm(temp_echelon_matrix[length(cols)+1,1:n]) < options.tol_nonzero_q    # so there is a circuit/stem to find in the current subset
                Z = nullspace(Vt[1:n,cols_plus])
                I = findall(x -> abs(x) > options.tol_coordinates, Z)                   # indices of the nonzero components - the real circuit
                stem = zeros(Int, size(Vt,2))
                value = dot(Z, Vt[n+1,cols_plus])                                       # quantity determining (a)symmetry
                if abs(value) > options.tol_nonzero_q                                   # stem of V, not Vt, so symmetric therefore requires the specific convention
                    stem[cols_plus[I]] = sign(value)*sign.(Z[I])                        # one has to take the right convention, given by this formula
                    append!(info.stems_sym_init, stem)                                  # update of info
                    info.nb_stems_sym += 1
                    # in this case, the recursion is continued but the info is modified (the recursion status changes)
                    # memory_bool becomes false, and the info of the linear combination found becomes the new one
                    isf_allsv_HnH_Wechelon_rec!(Vt, cols_plus, perm, temp_echelon_matrix, false, info, options)
                else                                                                    # actually a stem vector for Vt found directly, labelled as 'asym' + simple convention
                    stem[cols_plus[I]] = sign(Z[I[1]])*sign.(Z[I])                      # construction of the stem vector with convention
                    append!(info.stems_asym_init, stem)                                 # update of info
                    info.nb_stems_asym += 1
                end
            else                                                                        # continue the recursion, still independent even with the new one
                isf_allsv_HnH_Wechelon_rec!(Vt, cols_plus, perm, temp_echelon_matrix, memory, info, options)
            end

        else            # false means a linear relation has already been found, but not directly being a stem for Vt

            if norm(temp_echelon_matrix[length(cols):length(cols)+1, 1:n]) < 2*options.tol_nonzero_q
                # requires 2 linear dependence relations so last 2 lines checked the relations are combined directly; one could test if the second is fortunately perfect
                Z = nullspace(Vt[1:n, cols_plus])
                null_element = dot(Z[:,1], Vt[n+1,cols_plus]) * Z[:,2] - dot(Z[:,2], Vt[n+1,cols_plus]) * Z[:,1]
                I = findall(x -> abs(x) > options.tol_coordinates, null_element)        # indices of the nonzero components - the real circuit

                stem = zeros(Int, p)
                stem[cols_plus[I]] = sign(null_element[I[1]])*sign.(null_element[I])    # construction of the stem vector with convention
                append!(info.stems_asym_init, stem)                                     # update of info
                info.nb_stems_asym += 1
            else                                                                        # continue the recursion, still independent even with the new one
                isf_allsv_HnH_Wechelon_rec!(Vt, cols_plus, perm, temp_echelon_matrix, memory, info, options)
            end
        end
    end
end

# rational version 
"""
`isf_allsv_HnH_Wechelon_rec!(Vt, perm, cols, starting_echelon, info, options)`

Inspired from the computations of TOPCOM with the echelonned form. 

Rational version of `isf_allsv_from_indices_nH!`. Uses LinearAlgebraX. 
Recursive call on the subset size to possibly modify info. 
In each 'first' call, the size is 1 so it checks a subset of size 2. 
The echelon form is used purely on the submatrices, and nullspace is called
when the last line becomes zero. 
While the current recursion is in R^n, it is similar as in the nH case. 
When it is in R^{n+1}, the last line is already zero and the 
stem vector found wasn't one in R^{n+1} (otherwise it is already known). 
Thus when the _two_ last lines are zero, in R^n because this is where 
the echelon form is computed, necessarily in R^{n+1} there is a circuit,
which is obtained through the combination of both. 

The 'symmetric' (H) and 'asymmetric' (nH) versions are in their respective 
functions (files) `isf_allsv_H_Wechelon_rec!(.jl)`, `isf_allsv_nH_Wechelon_rec!(.jl)`.

For additional details, see `isf_allsv_from_indices_HnH(.jl)`, which is a more natural 
introduction, but often slightly slower. 
"""
function isf_allsv_HnH_Wechelon_rec_r!(Vt::Matrix, cols::Vector{Int64}, perm::Vector{Int64}, starting_echelon::Matrix, memory::Bool, info::Info, options::Options)

    p = size(Vt, 2)
    n = size(Vt, 1) - 1
    s = length(cols)

    for new_index in maximum(cols)+1:p                                                  # launch the following recursive calls
        cols_plus = [cols;new_index]
        temp_echelon_matrix = echelonned_form([starting_echelon ; Vt[1:n,new_index]'])  # echelonned matrix with the new vector

        if memory                                                                       # means previous vectors were still independent (natural phase, R^n)
            if norm(temp_echelon_matrix[s+1,1:n]) == 0                                  # so there is a circuit/stem to find in the current subset
                Z = nullspacex(Vt[1:n,cols_plus])                                       # use of LinearAlgebraX
                I = findall(x -> x != 0, Z)                                             # indices of the nonzero components - the real circuit
                stem = zeros(Int, size(Vt,2))
                value = dot(Z, Vt[n+1,cols_plus])                                       # quantity determining (a)symmetry
                if abs(value) > 0                                                       # stem of V, not Vt, so symmetric therefore requires the specific convention
                    stem[cols_plus[I]] = Int(sign(value))*Int.(sign.(Z[I]))             # the convention is different - one has to take the right convention, given by this formula
                    append!(info.stems_sym_init, stem)                                  # update of info
                    info.nb_stems_sym += 1
                    # in this case, the recursion is continued but the info is modified (the recursion status changes)
                    # memory_bool becomes false, and the info of the linear combination found becomes the new one
                    isf_allsv_HnH_Wechelon_rec_r!(Vt, cols_plus, perm, temp_echelon_matrix, false, info, options)
                else                                                                    # actually a stem vector for Vt found directly, labelled as 'asym' + simple convention
                    stem[cols_plus[I]] = Int(sign(Z[I[1]]))*Int.(sign.(Z[I]))           # construction of the stem vector with convention
                    append!(info.stems_asym_init, stem)                                 # update of info
                    info.nb_stems_asym += 1
                end
            else                                                                        # continue the recursion, still independent even with the new one
                isf_allsv_HnH_Wechelon_rec_r!(Vt, cols_plus, perm, temp_echelon_matrix, memory, info, options)
            end

        else            # false means a linear relation has already been found, but not directly being a stem for Vt

            if norm(temp_echelon_matrix[s:s+1, 1:n]) == 0       
                # requires 2 linear dependence relations so last 2 lines checked the relations are combined directly; one could test if the second is fortunately perfect
                Z = nullspacex(Vt[1:n, cols_plus])                                      # use of LinearAlgebraX
                null_element = dot(Z[:,1], Vt[n+1,cols_plus]) * Z[:,2] - dot(Z[:,2], Vt[n+1,cols_plus]) * Z[:,1]
                I = findall(x -> x != 0, null_element)                                  # indices of the nonzero components - the real circuit

                stem = zeros(Int, p)
                stem[cols_plus[I]] = Int(sign(null_element[I[1]]))*Int.(sign.(null_element[I])) # construction of the stem vector with convention
                append!(info.stems_asym_init, stem)                                     # update of info
                info.nb_stems_asym += 1

            else                                                                        # continue the recursion, still independent even with the new one
                isf_allsv_HnH_Wechelon_rec_r!(Vt, cols_plus, perm, temp_echelon_matrix, memory, info, options)
            end
        end
    end
end