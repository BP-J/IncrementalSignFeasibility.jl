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
`isf_sv_add_H!(element, perm, info, options)`

Add a stem vector to the list during heuristic D2. 
The convention is such that one has V[:,perm[1:nv]]*element = 0, 
where nv = length(element). 

The variants for cases nH and HnH are done separately (same file).
"""
function isf_sv_add_H!(element::Vector{Float64}, perm::Vector{Int64}, info::Info, options::Options)
    p = length(perm)

    I = findall(x -> abs(x) > options.tol_coordinates, element)             # indices of the nonzero components - the real circuit
    s = zeros(Int, p)
    s[perm[I]] = Int(sign(element[I[1]])) * Int.(sign.(element[I]))         # construction of the stem vector with convention

    # by construction, this stem vector is new - if it is not, then it has been found before, 
    # but then the second time the covering test should have found it
    info.stems_sym = [info.stems_sym s]                                     # update of info
    info.nb_stems_sym += 1
    append!(info.stem_sizes_sym, sum(abs.(s)))
end

#-----------------------------------------------------------------------

"""
`isf_sv_add_nH!(element, perm, type, info, options)`

Add a stem vector to the list during heuristic D2. 
The convention is such that one has V[:,perm[1:nv]]*element = 0, 
where nv = length(element). 

The variants for cases H and HnH are done separately (same file).
"""
function isf_sv_add_nH!(element::Vector{Float64}, type::String, perm::Vector{Int}, info::Info, options::Options)
    p = length(perm)

    I = findall(x -> abs(x) > options.tol_coordinates, element)             # indices of the nonzero components - the real circuit
    s = zeros(Int, p)

    # by construction, this stem vector is new - if it is not, then it has been found before, 
    # but then the second time the covering test should have found it
    if type == "sym"
        s[perm[I]] = Int(sign(element[I[1]])) * Int.(sign.(element[I]))     # construction of the stem vector with convention
        info.stems_sym = [info.stems_sym s]                                 # update of info
        info.nb_stems_sym += 1
        append!(info.stem_sizes_sym, sum(abs.(s)))
    else # type == "asym"
        s[perm[I]] = Int.(sign.(element[I]))                                # construction of the stem vector with convention
        info.stems_asym = [info.stems_asym s]                               # update of info
        info.nb_stems_asym += 1
        append!(info.stem_sizes_asym, sum(abs.(s)))
    end # end of the "sym" / "asym" disjonction
end

#-----------------------------------------------------------------------

"""
`isf_sv_add_HnH!(V, element, perm, type, info, options)`

Add a stem vector to the list during heuristic D2. 
The convention is such that one has V[:,perm[1:nv]]*element = 0, 
where nv = length(element). 

The variants for cases H and nH are done separately (same file).
"""
function isf_sv_add_HnH!(Vt::Matrix, element::Vector{Float64}, type::String, perm::Vector{Int}, info::Info, options::Options)
    p = length(perm)
    nv = length(element)
    
    I = findall(x -> abs(x) > options.tol_coordinates, element)                         # indices of the nonzero components - the real circuit
    s = zeros(Int, p)

    # by construction, this stem vector is new - if it is not, then it has been found before, 
    # but then the second time the covering test should have found it
    if type == "sym"
        s[perm[I]] = sign(dot(element, Vt[size(Vt,1), perm[1:nv]])) * sign.(element[I]) # construction of the stem vector with convention
        info.stems_sym = [info.stems_sym s]                                             # update of info
        info.nb_stems_sym += 1
        append!(info.stem_sizes_sym, sum(abs.(s)))
    else # type == "asym"
        s[perm[I]] = sign(element[I[1]]) * sign.(element[I])                            # the scalar product above is zero, so convention of the first nonzero sign
        info.stems_asym = [info.stems_asym s]                                           # update of info
        info.nb_stems_asym += 1
        append!(info.stem_sizes_asym, sum(abs.(s)))
    end # end of the "sym" / "asym" disjonction
end