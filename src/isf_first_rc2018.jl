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
`isf_first_rc2018!(Vt, info, options, values)`
Prepares the recursive call(s) to `isf_rec_rc2018`.
"""
function isf_first_rc2018!(Vt::Matrix, info::Info, options::Options, values::Values)

    p0 = size(Vt,2)
    Vt, colsel, colout = isf_noncolin!(Vt, options)
    # now the vectors & right hand sides are normed

    if options.symmetry
        n, p = size(Vt)
    else
        p = size(Vt, 2)
        n = size(Vt, 1) - 1
    end

    # other values are initialized by default to the right starting values

    info.flag = values.success

    if options.symmetry
        svec = [+1]
    else
        svecp = [+1]
        svecm = [-1]
    end

    perm = collect(1:p)
    # recursive run: the recursive process starts with svec == 1
    # in the symmetric case, the original code suffices
    # otherwise, at each node 

    if options.symmetry
        x = copy(Vt[1:n,1]) # the appropriate first vector
    else
        s1 = sign(Vt[n+1,1])
        if s1 > 0
            xp = 2 * Vt[n+1,1] * Vt[1:n,1] / (norm(Vt[1:n,1])^2)
            xm = zeros(n)
        elseif s1 < 0
            xp = zeros(n)
            xm = 2 * Vt[n+1,1] * Vt[1:n,1] / (norm(Vt[1:n,1])^2)
        else
            xp =  copy(Vt[1:n,1])
            xm = -copy(Vt[1:n,1])
        end
    end    

    # recursive call
    if p > 1
        if options.symmetry
            isf_rec_rc2018!(Vt, svec, perm, x, info, options, values)
        else
            isf_rec_rc2018!(Vt, svecp, perm, xp, info, options, values)
            isf_rec_rc2018!(Vt, svecm, perm, xm, info, options, values)
        end
        info.flag > 0 && return info
    end

#-------------------------------------------------------------------------------------------------------------------------------
# if some columns were colinear, they have been removed so they must be added again
#-------------------------------------------------------------------------------------------------------------------------------

    if p0 != p
        eliminated_columns = setdiff(1:p0, colsel)
        mat_sign_vec = zeros(Int, info.ns, p0)
        mat_sign_vec[:, colsel] = hcat(info.s...)'  # puts a list of vectors into a matrix (same 'dimensions')
        for elim in eliminated_columns
            if colout[elim] > 0
                mat_sign_vec[:,elim] = mat_sign_vec[:,abs(colout[elim])]
            else
                mat_sign_vec[:,elim] = -mat_sign_vec[:,abs(colout[elim])]
            end
        end
        info.s = [mat_sign_vec[i,:] for i in 1:info.ns]
    end
end

