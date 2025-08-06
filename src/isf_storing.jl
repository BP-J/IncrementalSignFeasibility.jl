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
`isf_storing!(svec, info, options, values)`

Storing of the current binary vector in the H / nH case
(the H case does not need the opposite and the nH does not know if it is valid).
"""
function isf_storing!(svec::Vector{Int}, p::Int, perm::Vector{Int}, info::Info)
    push!(info.s, zeros(Int, p))
    info.s[info.ns][perm] = svec
end

#-----------------------------------------------------------------------

"""
`isf_storing_HnH!(svec, HnH_info, symmetry_info, info, options, values)`

Storing of the current binary vector (and/or opposite) in the HnH case. 
"""
function isf_storing_HnH!(svec::Vector{Int}, p::Int, perm::Vector{Int}, symmetry_info::Int, info::Info)

    if symmetry_info == 0                   # the recurrence is still symmetric so both are added
        push!(info.s, zeros(Int,p))
        info.s[info.ns-1][perm] = svec 
        push!(info.s, zeros(Int,p))
        info.s[info.ns][perm] = - svec 
    else                                    # one adds only the correct one, using the symmetry info
        push!(info.s, zeros(Int, p))
        info.s[info.ns][perm] = symmetry_info*svec
    end
end
