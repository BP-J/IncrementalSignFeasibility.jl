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
`isf_rec_nH!(Vt, svec, perm, d, info, options, values)`

Recursive process forming the S-tree, using directions. 
Depending on the options, some improvements are added: 
QR factorization, bounds to avoid solving LO problems, 
heuristic on the next vector to add. 

The three versions are done separately, though the 
H and nH versions are very similar. 
"""
function isf_rec_nH!(Vt::Matrix, svec::Vector{Int}, perm::Vector{Int}, x::Vector, info::Info, options::Options, values::Values)

    info.flag = values.success
    p    = size(Vt,2)
    n    = size(Vt,1)-1
    nv   = length(svec)
    nvp  = nv+1

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
            I = findmax(abs.(Vt[1:n,Tc]' * x .- Vt[n+1,Tc]))[2]  # returns the first found
            ivp = Tc[I]                         # selected index
        elseif options.bestv == 2
            # makes no sense here
            println("no nH equivalent for bestv == 2, default -> perm[k+1]")
            ivp = perm[nvp]                         # selected index
        else            # == 3
            Tb = []     # list with a single descendent
            for i in Tc
                twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_nH(Vt, svec, x, Vt[:,i], T, options)
                if !twodesc
                    push!(Tb, i)
                end
            end
            if isempty(Tb)
                Tb = Tc
            end
            I = findmax(abs.(Vt[1:n,Tb]' * x .- Vt[n+1,Tb]))[2] # returns the first found
            ivp = Tb[I]
        end
        perm[nvp] = ivp         # memorize ivp in perm
    end

    vTx = Vt[1:n,ivp]' * x - Vt[n+1,ivp]

#-------------------------------------------------------------------------------------------------------------------------------
# Store the 'current' svsprod that will be modified for the recursive calls
#-------------------------------------------------------------------------------------------------------------------------------

    if options.svsprod
        svsprod_common_sym  = info.svsprod_sym
        svsprod_common_asym = info.svsprod_asym
    end

#-------------------------------------------------------------------------------------------------------------------------------
# Determine whether s has 1 or 2 descendents (depends on the options.dvnear0)
#-------------------------------------------------------------------------------------------------------------------------------

    if options.dvnear0 || abs(vTx) < options.tol_nonzero_q
        twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_nH(Vt, svec, x, Vt[:,ivp], T, options)
    else
        twodesc = false
    end

#-------------------------------------------------------------------------------------------------------------------------------
# vTx near zero ==> s has 2 descendents (without having to use an optimization problem to decide)
#-------------------------------------------------------------------------------------------------------------------------------

    if twodesc
        svec = [svec;+1]
        ## computations for the first one
        if minratiosJ == +Inf
            t = ratio+1
        else
            t = (ratio+minratiosJ)/2
        end
        xp = x + t*Vt[1:n,ivp]
        
        if nvp < p
            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym  = svsprod_common_sym  + info.stems_sym[ivp, :]
                info.svsprod_asym = svsprod_common_asym + info.stems_asym[ivp, :]
                info.cput_cover += tok()
            end
            isf_rec_nH!(Vt, svec, perm, xp, info, options, values)
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
        xm = x + t*Vt[1:n,ivp]

        if options.sv == 2 && options.svsprod # in this case, some stem vectors might be added so we (might) have to augment the current size
            if info.nb_stems_sym > length(svsprod_common_sym)   # stem added, so current size should be corrected
                @suppress begin
                    tick()
                end
                tempo = zeros(Int,p)
                tempo[perm[1:nv]] = svec[1:nv]
                svsprod_common_sym = [svsprod_common_sym ; info.stems_sym[:,length(svsprod_common_sym)+1:info.nb_stems_sym]'* tempo]
                info.cput_cover += tok()
            end
            if info.nb_stems_asym > length(svsprod_common_asym)   # stem added, so current size should be corrected
                @suppress begin
                    tick()
                end
                tempo = zeros(Int,p)
                tempo[perm[1:nv]] = svec[1:nv]
                svsprod_common_asym = [svsprod_common_asym ; info.stems_asym[:,length(svsprod_common_asym)+1:info.nb_stems_asym]' * tempo]
                info.cput_cover += tok()
            end
        end

        if nvp < p
            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym  = svsprod_common_sym  - info.stems_sym[ivp,:]
                info.svsprod_asym = svsprod_common_asym - info.stems_asym[ivp,:]
                info.cput_cover += tok()
            end
            isf_rec_nH!(Vt, svec, perm, xm, info, options, values)
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

#-------------------------------------------------------------------------------------------------------------------------------
# vTx != 0 ==> 1 single sign vector or 2 if optimization detects it
#-------------------------------------------------------------------------------------------------------------------------------

    else        # first sign is feasible, second checked
        svTx = Int(sign(vTx))     # the easy sign
        svec = [svec;svTx]

        if p > nvp 
            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym  += svTx*info.stems_sym[ivp,:]
                info.svsprod_asym += svTx*info.stems_asym[ivp,:]
                info.cput_cover += tok()
            end
            isf_rec_nH!(Vt, svec, perm, x, info, options, values)
            info.flag > 0 && return 
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

        # other sign checked
        svec[nvp] = -svTx

        if options.sv == 0
            covering_index = false         # no stem vectors so none of them is covered
        else
            if options.sv == 2 && options.svsprod  # in this case, some stem vectors might be added so we (might) have to augment the current size
                if info.nb_stems_sym > length(svsprod_common_sym)   # stem added, so current size should be corrected
                    @suppress begin
                        tick()
                    end
                    tempo = zeros(Int,p)
                    tempo[perm[1:nv]] = svec[1:nv]
                    svsprod_common_sym = [svsprod_common_sym ; info.stems_sym[:,length(svsprod_common_sym)+1:info.nb_stems_sym]'* tempo]
                    info.cput_cover += tok()
                end
                if info.nb_stems_asym > length(svsprod_common_asym)   # stem added, so current size should be corrected
                    @suppress begin
                        tick()
                    end
                    tempo = zeros(Int,p)
                    tempo[perm[1:nv]] = svec[1:nv]
                    svsprod_common_asym = [svsprod_common_asym ; info.stems_asym[:,length(svsprod_common_asym)+1:info.nb_stems_asym]' * tempo]
                    info.cput_cover += tok()
                end
            end
            info.nb_sv_checks += 1

            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym  = svsprod_common_sym  - svTx*info.stems_sym[ivp,:]
                info.svsprod_asym = svsprod_common_asym - svTx*info.stems_asym[ivp,:]
                covering_index = isf_cover_stem_nH_bis(svec, info, options)
                info.cput_cover = tok()
            else
                @suppress begin
                    tick()
                end
                covering_index = isf_cover_stem_nH(svec, perm, info, options)
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
            x, lambda = isf_feas!(svec[1:nv]' .* Vt[:,T], -svTx*Vt[1:n,ivp], -svTx*Vt[n+1,ivp], info, options, values)
            info.cput_lp += tok()
            info.flag > 0 && return 

            if any(x .!= Inf)      # system is feasible
                if p > nvp
                    isf_rec_nH!(Vt, svec, perm, x, info, options, values)
                    info.flag > 0 && return 
                else
                    info.ns += 1
                    options.s && isf_storing!(svec, p, perm, info)
                end
            else
                if options.sv == 2      # dual variable becomes a stem vector
                    if abs((svec.* lambda[1:nvp])' * Vt[n+1,[T;ivp]]) < options.tol_nonzero_q
                        isf_sv_add_nH!(svec.* lambda[1:nvp], "sym", perm, info, options)
                    else
                        isf_sv_add_nH!(svec.* lambda[1:nvp], "asym", perm, info, options)
                    end
                end
            end
        end
    end     # of the vTx != 0 test
end