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

### The instances generated here are taken from 
### "COMPUTING CHARACTERISTIC POLYNOMIALS OF HYPERPLANE ARRANGEMENTS WITH SYMMETRIES", by
### TAYLOR BRYSIEWICZ, HOLGER EBLE, AND LUKAS KÜHNE
### https://arxiv.org/pdf/2105.14542.pdf

"""
`matrix = isf_generation_threshold(d)`
Returns the matrix corresponding to the threshold arrangement in dimension d(+1).
The 2^d vectors are the vectors (1, c_1, ..., c_d) for c_i ∈ {0,1}, computed recursively.
"""
function isf_generation_threshold(d)
    M = [0 1]
    for i in 2:d
        M = [[M ; zeros(Int, 1,2^(i-1))] [M ; ones(Int, 1,2^(i-1))]]
    end
    return [ones(Int, 1, 2^d) ; M]
end

"""
`matrix = isf_generation_resonance(d)`
Returns the matrix correspondign to the threshold arrangement in dimension d(+1).
The 2^d - 1 vectors are the vectors (1, c_1, ..., c_d) for c_i ∈ {0,1} and c != 0_d.
Uses the threshold arrangement.
"""
function isf_generation_resonance(d)
    M = isf_generation_threshold(d)
    return M[2:d+1, 2:2^d]
end

"""
`matrix = isf_generation_crosspoly(d)`
Returns the matrix corresponding to the cross-polytope arrangement in dimension d(+1).
The 2d vectors are the (1, ± e_i) for e_i the i-th canonical vector.
"""
function isf_generation_crosspoly(d)
    M = ones(Int, 1, 2*d)
    M = [M ; I -I]
    return M
end

### the permutohedron ones are different from the RC permutahedron instances, therefor not done. 

"""
`matrix = isf_generation_demicube(d)`
Rreturns the matrix corresponding to the demicube arrangement in dimension d(+1).
The 2^(d-1) vectors are the subsets of the threshold arrangements having odd (even with the additional 1) number of 1's.
Done with the powerset function.
"""
function isf_generation_demicube(d)
    M = zeros(Int, d, 2^(d-1))
    base = [i for i in 1:d]
    count = 0
    for i in 1:2:d
        indices = collect( powerset(base, i, i) )

        for j in eachindex(indices)
            M[indices[j], count+j] .+= 1
        end
        count += length(indices)
    end
    return [ones(Int, 1, 2^(d-1)) ; M]
end
