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
# Computations related to modification B.
# Nodes of the S-tree always have one descendant, sometimes two
# answering this question is the main goal, but instead of using 
# a linear optimization problem sometimes simple computations can yield the result

# The computations are rather simple, short and similar for all
# versions therefore they are grouped in this common file. 

"""
`twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_H(V, s, d, v, T, options)`

Returns true if there are two easy descendants, as well as the associated 
relevant values. Code associated to the modification B, i.e. elaborated 
bounds/tests to get two descendants without LO solve.
Version for symmetric arrangements. 
"""
function isf_ntype_H(V::Matrix, s::Vector{Int64}, d::Vector, v::Vector, T::Vector{Int64}, options::Options)

    # quantities s_i (v_i, v) to know whether the indices belong in the left or right set
    sVTv = s .* (V[:,T]' * v)

    # verifying if there are 2 descendants by computing the quantities
    # -(v_i, d)/ (v, v_i)
    ratio  = -v'*d / (v'*v)
    ratios = (-V[:,T]' * d) ./ (V[:,T]' * v)
    I = findall(y -> y > 0, sVTv)
    J = findall(y -> y < 0, sVTv)

    # comparing with the indices such that s_i (v_i, v) > 0
    if isempty(I)
        maxratiosI = -Inf
        testI = true
    else
        maxratiosI = maximum(ratios[I])
        testI = (maxratiosI < ratio - options.tol_nonzero_q)
    end

    # comparing with the indices such that s_i (v_i, v) < 0
    if isempty(J)
        minratiosJ = +Inf
        testJ = true
    else
        minratiosJ = minimum(ratios[J])
        testJ = (ratio + options.tol_nonzero_q < minratiosJ)
    end
    return testI && testJ, ratio, maxratiosI, minratiosJ
end

#-----------------------------------------------------------------------

"""
`twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_nH(V, s, x, v, T, options)`

Returns true if there are two easy descendants, as well as the associated 
relevant values. Code associated to the modification B, i.e. elaborated 
bounds/tests to get two descendants without LO solve.
Version for arbitrary arrangements. 
"""
function isf_ntype_nH(Vt::Matrix, s::Vector{Int64}, x::Vector, v::Vector, T::Vector{Int64}, options::Options)

    n = size(Vt,1)-1

    # quantities s_i (v_i, v) to know whether the indices belong in the left or right set
    sVTv = s .* (Vt[1:n,T]'* v[1:n])

    # verifying if there are 2 descendants by computing the quantities
    # [\tau_i -(v_i, x)]/ (v, v_i)
    ratios = (Vt[n+1,T] .- Vt[1:n,T]' * x) ./ (Vt[1:n,T]' * v[1:n])
    ratio  = (v[n+1] - x'*v[1:n])/(v[1:n]' * v[1:n])

    I = findall(y -> y > 0, sVTv)
    J = findall(y -> y < 0, sVTv)

    # comparing with the indices such that s_i (v_i, v) > 0
    if isempty(I)
        maxratiosI = -Inf
        testI = true
    else
        maxratiosI = maximum(ratios[I])
        testI = (maxratiosI < ratio - options.tol_nonzero_q)
    end

    # comparing with the indices such that s_i (v_i, v) < 0
    if isempty(J)
        minratiosJ = +Inf
        testJ = true
    else
        minratiosJ = minimum(ratios[J])
        testJ = (ratio + options.tol_nonzero_q < minratiosJ)
    end
    return testI && testJ, ratio, maxratiosI, minratiosJ
end

#-----------------------------------------------------------------------

"""
`twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_HnH(V, s, d, v, T, info, options)`

Returns true if there are two easy descendants, as well as the associated 
relevant values. Code associated to the modification B, i.e. elaborated 
bounds/tests to get two descendants without LO solve. The equations/formulae 
are slightly modified in the nH case. Version for the compact algorithms.
The parameter info describes the status of the HnH recursion. 
"""
function isf_ntype_HnH(Vt::Matrix, s::Vector{Int64}, d::Vector, v::Vector, T::Vector{Int64}, info::Int, options::Options)

    # the same computations are done but the size is adapted
    if info == 0 # 
        N = size(Vt,1)-1
    else
        N = size(Vt,1)
    end
    v = v[1:N]
    Vt = Vt[1:N,:]

    # quantities s_i (v_i, v) to know whether the indices belong in the left or right set
    sVTv = s .* (Vt[:,T]' * v)

    # verifying if there are 2 descendants by computing the quantities
    # -(v_i, d)/ (v, v_i)
    ratio  = -v'*d / (v'*v)
    ratios = (-Vt[:,T]' * d) ./ (Vt[:,T]' * v)
    I = findall(y -> y > 0, sVTv)
    J = findall(y -> y < 0, sVTv)

    # comparing with the indices such that s_i (v_i, v) > 0
    if isempty(I)
        maxratiosI = -Inf
        testI = true
    else
        maxratiosI = maximum(ratios[I])
        testI = (maxratiosI < ratio - options.tol_nonzero_q)
    end

    # comparing with the indices such that s_i (v_i, v) < 0
    if isempty(J)
        minratiosJ = +Inf
        testJ = true
    else
        minratiosJ = minimum(ratios[J])
        testJ = (ratio + options.tol_nonzero_q < minratiosJ)
    end
    return testI && testJ, ratio, maxratiosI, minratiosJ
end