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
`info = isf_rec_rc2018!(V, svec, perm, d, info, options, values)`

Recursive procedure, forming the S-tree.
"""
function isf_rec_rc2018!(V::Matrix, svec::Vector{Int}, perm::Vector{Int64}, x::Vector, info::Info, options::Options, values::Values)

    info.flag = values.success

    p = size(V, 2)      # total number of vectors
    nv = length(svec)   # number of vectors associated with svec
    nvp = nv+1

    if options.symmetry
        n = size(V,1)   # either symmetric and told by user or symmetric and useless last line of zeros removed
        rhs = 0         # for optimization purposes below
    else
        n = size(V,1) - 1
        rhs = V[n+1,nvp]
    end

    v   = V[1:n,nvp]
    vTx = v' * x

    if vTx - rhs > options.tol_nonzero_q
        svec = [svec;+1]

        if p > nvp 
            isf_rec_rc2018!(V, svec, perm, x, info, options, values)
            info.flag > 0 && return info
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

        # LO solved to determine if the other descendant

        svec[nvp] = -1

        @suppress begin
            time_1 = tick()
        end

        x, lambda = isf_feas!(svec[1:nv]' .* V[:,1:nv], -v, -rhs, info, options, values)
        
        info.cput_lp += tok()

        info.flag > 0 && return info
        if any(x .!= Inf) # feasible system, if infeasible Inf * ones(n) is returned
            if p > nvp 
                isf_rec_rc2018!(V,svec,perm,x,info,options,values)
                info.flag > 0 && return info
            else
                info.ns += 1
                options.s && isf_storing!(svec, p, perm, info)
            end
        end

    elseif vTx - rhs < -options.tol_nonzero_q # on the other side of the hyperplane, so the other sense
        svec = [svec;-1]

        if p > nvp 
            isf_rec_rc2018!(V, svec, perm, x, info, options, values)
            info.flag > 0 && return info
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

        # LO solved to determine if the other descendant

        svec[nvp] = +1

        @suppress begin
            time_1 = tick()
        end
        
        x, lambda = isf_feas!(svec[1:nv]' .* V[:,1:nv], +v, +rhs, info, options, values)
        
        info.cput_lp += tok()

        info.flag > 0 && return info
        if any(x .!= Inf) # feasible system
            
            if p > nvp 
                isf_rec_rc2018!(V,svec,perm,x,info,options,values)
                info.flag > 0 && return info
            else
                info.ns += 1
                options.s && isf_storing!(svec, p, perm, info)
            end
        end
    else # vTx == rhs, the annoying case (the == 0 is only in this version, other have an appropriate test) : both subtrees are continued and a LOP is artificially added 

        if p > nvp
            sVTv = svec .* (V[1:n,1:nv]' * v)
            if options.symmetry
                ratios = (-V[:,1:nv]' * x) ./ (V[:,1:nv]' * v)
            else
                ratios = (V[n+1,1:nv] - V[1:n,1:nv]' * x) ./ (V[1:n,1:nv]' * v)
            end

            I = findall(y -> y > 0, sVTv)
            J = findall(y -> y < 0, sVTv)
            # left part
            if isempty(I)
                maxratiosI = -2
            else
                maxratiosI = maximum(ratios[I])
            end

            # right part
            if isempty(J)
                minratiosJ = +2
            else
                minratiosJ = minimum(ratios[J])
            end
            isf_rec_rc2018!(V, [svec;+1], perm, x + 0.5*minratiosJ*v, info, options, values)
            info.flag > 0 && return info
            isf_rec_rc2018!(V, [svec;-1], perm, x + 0.5*maxratiosI*v, info, options, values)
            info.flag > 0 && return info

        else
            info.ns += 2
            options.s && (info.ns -= 1; isf_storing!([svec;+1], p, perm, info); info.ns += 1; isf_storing!([svec;-1], p, perm, info))
        end

        # no need to use the feasibility part, done to artificially correspond to original RC
        info.nb_losolve += 1
        info.nb_feaslop += 1
    end
end
