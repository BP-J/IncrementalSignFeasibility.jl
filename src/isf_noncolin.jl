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
`W, colsel, colout = isf_noncolin!(V, info, options, values)`

Removes nonlinear columns and norms the matrix.
Selects the noncolinear columns of V (Vt) ones, assuming they are normalized.
The second output represents the selected columns, while the third
is an array for which [j] is the index in colsel of V[:,j] (j <= p).
This function need not change for the nonhomogeneous versions, since
hyperplanes i and j are identical if V[:,i] and V[:,j] (Vt[:,i] and Vt[:,j])
are colinear in homogeneous (nonhomogeneous) instances. 
"""
function isf_noncolin!(V::Matrix, options::Options)
    p = size(V,2)

    norms = sqrt.(sum(V .* V, dims=1))  # norming the data
    V = V ./ norms

    colsel = collect(1:p)               # initially all columns are supposed noncolinear
    colout = copy(colsel)               # pointers to colsel telling where are the columns in the future W (> 0 if same 'sense')

    js = 0
    for jo in 1:p-1
        if colout[jo]< jo               # column jo has already been detected as colinear to another
            continue
        end
        js += 1                         # one more noncolinear column
        colsel[js] = jo                 # selected in V, index jo
        colout[jo] = js                 # the column 'jo' of V is the 'js' one in colsel
        v = V[:,jo]
        for j = jo+1:p
            if colout[j] < j            # skipped if coulout[j] has already been seen as colinear
                continue
            end                         # since normed, vectors are colinear if they are equal or opposite (up to tolerance)
            if norm(v - V[:,j], Inf) <= options.tol_coordinates
                colout[j] = jo
            elseif norm(v + V[:,j], Inf) <= options.tol_coordinates
                colout[j] = -jo
            end
        end
    end

    if colout[p] == p                   # the last one is not colinear to another
        js = js+1
        colsel[js] = p
        colout[p] = js
    end

    W = V[:, colsel[1:js]]              # the matrix used in the rest
    colsel = colsel[1:js]

    return W, colsel, colout
end

# rational version
"""
`W, colsel, colout = isf_noncolin_r!(V, info, options, values)`

Removes nonlinear columns and norms the matrix, rational version. 
Selects the noncolinear columns of V (Vt) ones, assuming they are normalized.
The second output represents the selected columns, while the third
is an array for which [j] is the index in colsel of V[:,j] (j <= p).
This function need not change for the nonhomogeneous versions, since
hyperplanes i and j are identical if V[:,i] and V[:,j] (Vt[:,i] and Vt[:,j])
are colinear in homogeneous (nonhomogeneous) instances. 
"""
function isf_noncolin_r!(V::Matrix)
    p = size(V,2)

    norms = maximum(abs, V, dims=1)     # norming the data with the infinite norm (one cannot use the 2-norm since its not rational, but one could use the 1-norm)
    V = V ./ norms

    colsel = collect(1:p)               # initially all columns are supposed noncolinear
    colout = copy(colsel)               # pointers to colsel telling where are the columns in the future W (> 0 if same 'sense')

    js = 0
    for jo in 1:p-1
        if colout[jo]< jo               # column jo has already been detected as colinear to another
            continue
        end
        js += 1                         # one more noncolinear column
        colsel[js] = jo                 # selected in V, index jo
        colout[jo] = js                 # the column 'jo' of V is the 'js' one in colsel
        v = V[:,jo]
        for j = jo+1:p
            if colout[j] < j            # skipped if coulout[j] has already been seen as colinear
                continue
            end
            if v == V[:,j]              # since normed and rational, two vectors are colinear if they are equal or opposite
                colout[j] = jo
            elseif v == -V[:,j]
                colout[j] = -jo
            end
        end
    end

    if colout[p] == p                   # the last one is not colinear to another
        js = js+1
        colsel[js] = p
        colout[p] = js
    end

    W = V[:, colsel[1:js]]              # the matrix used in the rest
    colsel = colsel[1:js]

    return W, colsel, colout
end