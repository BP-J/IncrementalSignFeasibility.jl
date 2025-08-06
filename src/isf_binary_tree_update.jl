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
# isf_binary_tree_update(binary_vec, pointer_vec)

# Shenanigans on binary vectors: manual update to the 'next' binary values
# according to the order used.

"""
`b, p = isf_binary_tree_update(b, p)`
b is the current binary vector (~ subset) and p the index pointer.
Moves to the next element to be studied in the binary tree.
Can also stop the exploration.
"""
function isf_binary_tree_update(b::Vector{Int64}, p::Int64)

    pmax = length(b)
    while true
        if p < pmax
            p += 1
            b[p] = 1
            return b, p
        end

        while (p > 0) && (b[p] == 0)    # backtrack to find the first 1
            p = p-1
        end
        if p == 0                       # end of the exploration
            return b, p
        end
        b[p] = 0                        # replaced by 0 and the main loop continues
    end
    return b, p
end