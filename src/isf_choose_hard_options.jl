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
`options = isf_choose_wechelon_option(V, Q, R, info, options)`

Empirical rule modifying the options.wechelon value.
One wants options.wechelon to be true if the circuits are large enough
in size so that building recursively the echelonned matrix is relevant. 
This is mainly done by looking at the QR factorization (not only the rank), 
checking if some parts of the matrix are (close to) zero. 
(The rank by itself is not extremely useful to describe the matrix.)
"""

function isf_choose_wechelon_option!(V::Matrix, Q, R, info::Info, options::Options)
    
    p = size(V,2)
    r = info.r

    if r <= 3                                           # the rank is small enough so the circuits will be too
        options.wechelon = false                        # therefore no recursive computation
    else
        n = size(Q,1)                               
        l = []
        for i in 1:n
            if norm(R[i,r+1:p]) < options.tol_nonzero_q # if true, then a part of R is zero, which, by the QR factorization, suggests there might be a subarrangement 
                l = append!(l,i)                        # of small intrinsic dimension
            end
        end
        rest = setdiff(1:n,l)
        c = 0
        for i in 1:n                                    
            if norm(Q[i,rest]) < options.tol_nonzero_q  # then one checks if Q has zeros on the complementary indices rest = [1:n] \ l
                c += 1                                  # c counts the number of zero lines in V[:,r+1:p]
            end
        end

        if r-c <= 3                                     # finally, if r-c <= 3, then there is a subarrangement in dimension <= 3, so no recursive computations
            options.wechelon = false
        else
            options.wechelon = true
        end
    end
end

# it might be optimal to do recursive for one stem type (asymmetric) and regular for the other (symmetric), but this has not been implemented. 

"""
`options = isf_choose_svsprods_option!(V, info, options)`

Empirical rule modifying the options.svsprod value.
For options.svsprod, first for D1 and D2 it is (quite naturally) bad.
However, for D3 and D4, the question is relevant. 
Observations on the tested instances suggested that for low amounts of stem vectors, 
therefore low rank and/or number of vectors, it is (logically) not worth.
A proposed test is to compare the number of stem vectors to the maximal number (general position), 
and thus to judge whether to use the recursive computation or not. 
"""
function isf_choose_svsprods_option!(V::Matrix, info::Info, options::Options)
    
    p = size(V,2)
    r = info.r
    
    if options.sv < 3                                           # not all stem vectors so no recursive computation
        options.svsprod = false
    else
        if r <= 3 || p <= 15                                    # if rank or number of vectors is small, no recursive computation ("few" stem vectors)
            options.svsprod = false
        elseif p >= 25                                          # assumed the number of stem vectors will be large anyways
            options.svsprod = true
        else
            if options.algorithm <= 7                           # split on whether it is an HnH algorithm or not
                if 2*(info.nb_stems_sym + info.nb_stems_asym) > binomial(p, r+1)
                    options.svsprod = true                      # assumed there is a large number of stem vectors, so recursive computation 
                else
                    options.svsprod = false
                end
            else                                                # HnH version, so a lot of both types of stem vectors
                if (2*info.nb_stems_sym > binomial(p, r+1)) || (2*info.nb_stems_asym > binomial(p, r+2))
                    options.svsprod = true                      # assumed there is a large number of stem vectors, so recursive computation 
                else
                    options.svsprod = false
                end
            end
        end
    end
end

