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
# isf_bin_operations contains very simple functions that update binary vectors. 
# isf_bin_minus is used in the H part of the code, 
# while for convenience purposes the nH code uses the isf_bin_plus function. 

"""
`b = isf_binary_minus(b)`
Binary decrease: makes the zeros to 1 and the first 1 to zero.
"""
function isf_binary_minus(b::Vector{Int64})

    if all(b .== 0) 
        return []
    end
    n = length(b)
    for i in 1:n
        if b[i] == 0
            b[i] = 1
        else
            b[i] = 0
            break
        end
    end
    return b
end

"""
`b = isf_binary_plus(b)`
Binary increase: makes the ones to zero and the first zero to 1.
"""
function isf_binary_plus(b::Vector{Int64})

    if all(b .== 1)
        return []
    end
    n = length(b)
    for i in 1:n
        if b[i] == 1
            b[i] = 0
        else
            b[i] = 1
            break
        end
    end
    return b
end

### Versions manipulating sign vectors instead of binary vectors. 

"""
`b = isf_decrement_s(s)`
Binary decrease: makes the zeros to 1 and the first 1 to zero
for a sign vector.
"""
function isf_decrement_s(s::Vector{Int64})

    if all(s .== -1) 
        return []
    end

    i = findfirst(t -> t == 1, s)
    s[i] = -1
    s[1:i-1] .= 1
    return s
end

"""
`b = isf_increment_s(s)`
Binary increase: makes the ones to zero and the first zero to 1
for a sign vector.
"""
function isf_increment_s(s::Vector{Int64})

    if all(s .== 1)
        return []
    end

    i = findfirst(t -> t == -1, s)
    s[i] = 1
    s[1:i-1] .= -1
    return s
end