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
`isf_rec_H!(V, svec, perm, d, info, options, values)`

Recursive process forming the S-tree, using directions. 
Depending on the options, some improvements are added: 
QR factorization, bounds to avoid solving LO problems, 
heuristic on the next vector to add. 

The three versions are done separately, though the 
H and nH versions are very similar. 
"""
function isf_rec_H!(V::Matrix, svec::Vector{Int}, perm::Vector{Int}, d::Vector, info::Info, options::Options, values::Values)
    
    info.flag = values.success
    p   = size(V,2)
    nv  = length(svec)
    nvp = nv+1

#-------------------------------------------------------------------------------------------------------------------------------
# Determine the next v; this operation must be done for each {s_i v_i: i in T}, since although the signs s_i play no role (with
# the adopted technique), d may change (and this has an impact on the choice of the new vector with the adopted techniques)
#-------------------------------------------------------------------------------------------------------------------------------

    T = perm[1:nv]          # index set of vectors used so far

    # 'ivp' will be the next vector to add

    if options.bestv == 0 
        ivp = perm[nvp]     # following the permutation from the QR factorization
    else
        Tc = setdiff(collect(1:p), T)
        if options.bestv == 1
            I = findmax(abs.(V[:,Tc]' * d))[2]  # returns the first found
            ivp = Tc[I]                         # selected index
        elseif options.bestv == 2
            dd = sum(V[:,T], dims=2)
            I = findmax(abs.(V[:,Tc]' * dd))[2] # returns the first found
            ivp = Tc[I]                         # selected index
        else            # == 3
            Tb = []     # list with a single descendent
            for i in Tc
                twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_H(V, svec, d, V[:,i], T, options)
                if !twodesc
                    push!(Tb, i)
                end
            end
            if isempty(Tb)
                Tb = Tc
            end
            I = findmax(abs.(V[:,Tb]' * d))[2] # returns the first found
            ivp = Tb[I]
        end
        perm[nvp] = ivp         # memorize ivp in perm
    end

    vTd = V[:,ivp]' * d

#-------------------------------------------------------------------------------------------------------------------------------
# Store the 'current' svsprod that will be modified for the recursive calls
#-------------------------------------------------------------------------------------------------------------------------------

    if options.svsprod
        svsprod_common = info.svsprod_sym
    end
    
#-------------------------------------------------------------------------------------------------------------------------------
# Determine whether s has 1 or 2 descendents (depends on the options.dvnear0)
#-------------------------------------------------------------------------------------------------------------------------------

    if options.dvnear0 || abs(vTd) < options.tol_nonzero_q
        twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_H(V, svec, d, V[:,ivp], T, options)
    else
        twodesc = false
    end

#-------------------------------------------------------------------------------------------------------------------------------
# vTd near zero ==> s has 2 descendents (without having to use an optimization problem to decide)
#-------------------------------------------------------------------------------------------------------------------------------
    
    if twodesc
        svec = [svec;+1]
        ## computations for the first one
        if minratiosJ == +Inf
            t = ratio+1
        else
            t = (ratio+minratiosJ)/2
        end
        dp = d + t*V[:,ivp]
        dp = dp/norm(dp)

        if p > nvp     
            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym = svsprod_common + info.stems_sym[ivp,:]
                info.cput_cover += tok()
            end

            isf_rec_H!(V, svec, perm, dp, info, options, values)
            info.flag > 0 && return
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

        svec[nvp] = -1
        ## computations for the second one
        if maxratiosI == -Inf
            t = ratio-1
        else
            t = (ratio+maxratiosI)/2
        end
        dm = d + t*V[:,ivp]
        dm = dm/norm(dm)

        if options.sv == 2 && options.svsprod # in this case, some stem vectors might be added so we (might) have to augment the current size
            if info.nb_stems_sym > length(svsprod_common)   # stem added, so current size should be corrected
                @suppress begin
                    tick()
                end
                tempo = zeros(Int,p)
                tempo[perm[1:nv]] = svec[1:nv]
                svsprod_common = [svsprod_common ; info.stems_sym[:,length(svsprod_common)+1:info.nb_stems_sym]'* tempo]
                info.cput_cover += tok()
            end
        end
        if p > nvp
            if options.svsprod      # update of the recursive computation, here as info.svsprod was modified, one must use the 'common' info stored at the beginning
                @suppress begin
                    tick()
                end
                info.svsprod_sym = svsprod_common - info.stems_sym[perm[nvp],:]
                info.cput_cover += tok()
            end

            isf_rec_H!(V, svec, perm, dm, info, options, values)
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

#-------------------------------------------------------------------------------------------------------------------------------
# vTd != 0 ==> 1 single sign vector or 2 if optimization detects it
#-------------------------------------------------------------------------------------------------------------------------------

    else        # first sign is feasible, second checked
        svTd = Int(sign(vTd))     # the easy sign
        svec = [svec;svTd]

        if p > nvp 
            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym += svTd*info.stems_sym[ivp,:]
                info.cput_cover += tok()
            end

            isf_rec_H!(V, svec, perm, d, info, options, values)
            info.flag > 0 && return
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end
        

        # other sign checked 
        svec[nvp] = -svTd

        if options.sv == 0
            covering_index = false         # no stem vectors so none of them is covered
        else
            if options.sv == 2 && options.svsprod # in this case, some stem vectors might be added so we (might) have to augment the current size
                if info.nb_stems_sym > length(svsprod_common)   # stem added, so current size should be corrected
                    @suppress begin
                        tick()
                    end
                    tempo = zeros(Int,p)
                    tempo[perm[1:nv]] = svec[1:nv]
                    svsprod_common = [svsprod_common ; info.stems_sym[:,length(svsprod_common)+1:info.nb_stems_sym]'* tempo]
                    info.cput_cover += tok()
                end
            end
            info.nb_sv_checks += 1

            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym = svsprod_common - svTd*info.stems_sym[ivp,:]
                covering_index = isf_cover_stem_H_bis(svec, info, options)
                info.cput_cover += tok()
            else
                @suppress begin
                    tick()
                end
                covering_index = isf_cover_stem_H(svec, perm, info, options)
                info.cput_cover += tok()
            end
        end

        if covering_index
            ## this means the current svec is infeasible
            info.nb_sv_detect += 1

        else        # a LO problem is solved to check [svec;-1]
            @suppress begin
                tick()
            end
            d, lambda = isf_feas!(svec[1:nv]' .* V[:,T], -svTd*V[:,ivp], 0, info, options, values)
            info.cput_lp += tok()
            info.flag > 0 && return

            if any(d .!= Inf)      # system is feasible
                if p > nvp
                    isf_rec_H!(V, svec, perm, d, info, options, values)
                    info.flag > 0 && return
                else
                    info.ns += 1
                    options.s && isf_storing!(svec, p, perm, info)
                end
            else
                if options.sv == 2      # dual variable becomes a stem vector
                    isf_sv_add_H!(svec.* [lambda[1:nv];1], perm, info, options)
                end
            end
        end
    end     # of the vTd != 0 test
end