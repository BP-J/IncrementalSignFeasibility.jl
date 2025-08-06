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
`isf_rec_HnH!(V, svec, perm, d, info, options, values)`

Recursive process forming the S-tree, using directions. 
Depending on the options, some improvements are added: 
QR factorization, bounds to avoid solving LO problems, 
heuristic on the next vector to add. 

The three versions are done separately, especially
because the HnH version is more complicated and uses 
additional inputs HnH_info, symmetry_info, which tell 
the current state of the recursion. 
- HnH_info == 0 means the 'natural' state, (R^n)
- HnH_info == 1 means the 'upgraded' state, i.e.,
  the current branch is V-infeasible but Vt-feasible
  symmetry_info is a sign telling which sign vector is really in S
(one could use only symmetry_info and "if symmetry_info != 0"...)
"""
function isf_rec_HnH!(Vt::Matrix, svec::Vector{Int}, perm::Vector{Int}, d::Vector, HnH_info::Int, symmetry_info::Int, info::Info, options::Options, values::Values)

    info.flag = values.success
    p    = size(Vt,2)
    n    = size(Vt,1) - 1 + HnH_info      # if HnH_info == 1, then size is n, 
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
            I = findmax(abs.(Vt[1:n,Tc]' * d))[2]  # returns the first found
            ivp = Tc[I]                         # selected index
        elseif options.bestv == 2
            dd = sum(V[1:n,T], dims=2)
            I = findmax(abs.(V[:,Tc]' * dd))[2] # returns the first found
            ivp = Tc[I]                         # selected index
        else            # == 3
            Tb = []     # list with a single descendent
            for i in Tc
                twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_HnH(Vt, svec, d, Vt[:,i], T, HnH_info, options)
                if !twodesc
                    push!(Tb, i)
                end
            end
            if isempty(Tb)
                Tb = Tc
            end
            I = findmax(abs.(Vt[1:n,Tb]' * d))[2] # returns the first found
            ivp = Tb[I]
        end

        perm[nvp] = ivp         # memorize ivp in perm

    end

    vTd = Vt[1:n,ivp]' * d

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

    if options.dvnear0 || abs(vTd) < options.tol_nonzero_q
        twodesc, ratio, maxratiosI, minratiosJ = isf_ntype_HnH(Vt, svec, d, Vt[:,ivp], T, HnH_info, options)
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
        dp = d + t*Vt[1:n,ivp]
        dp = dp/norm(dp)
        
        if nvp < p
            if options.svsprod      # update of the recursive computation
                @suppress begin
                    tick()
                end
                info.svsprod_sym  = svsprod_common_sym  + info.stems_sym[ivp, :]
                info.svsprod_asym = svsprod_common_asym + info.stems_asym[ivp, :]
                info.cput_cover += tok()
            end
            isf_rec_HnH!(Vt, svec, perm, dp, HnH_info, symmetry_info, info, options, values)
            info.flag > 0 && return 
        else
            info.ns += 2-HnH_info
            options.s && isf_storing_HnH!(svec, p, perm, symmetry_info, info)
        end

        svec[nvp] = -1
        ## computations for the second one
        if maxratiosI == -Inf
            t = ratio-1
        else
            t = (ratio+maxratiosI)/2
        end
        dm = d + t*Vt[1:n,ivp]

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
            isf_rec_HnH!(Vt, svec, perm, dm, HnH_info, symmetry_info, info, options, values)
        else
            info.ns += 2-HnH_info
            options.s && isf_storing_HnH!(svec, p, perm, symmetry_info, info)
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
                info.svsprod_sym  += svTd*info.stems_sym[ivp,:]
                info.svsprod_asym += svTd*info.stems_asym[ivp,:]
                info.cput_cover += tok()
            end
            isf_rec_HnH!(Vt, svec, perm, d, HnH_info, symmetry_info, info, options, values)
            info.flag > 0 && return 
        else
            info.ns += 2-HnH_info
            options.s && isf_storing_HnH!(svec, p, perm, symmetry_info, info)
        end

        # other sign checked
        svec[nvp] = -svTd

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
                info.svsprod_sym  = svsprod_common_sym  - svTd*info.stems_sym[ivp,:]
                info.svsprod_asym = svsprod_common_asym - svTd*info.stems_asym[ivp,:]
                covering_index, detail, sgn = isf_cover_stem_HnH_bis(svec, HnH_info, info, options)
                info.cput_cover += tok()
            else
                @suppress begin
                    tick()
                end
                covering_index, detail, sgn = isf_cover_stem_HnH(svec, perm, HnH_info, info, options)
                info.cput_cover += tok()
            end
        end

        if covering_index
            ## this means the current svec is infeasible
            info.nb_sv_detect += 1
            if detail == "asym"     # means the current sign vector is infeasible even in dimension n+1, so the recursion is stopped

            else        # covering_index != [] means detail ∈ ["sym", "asym"], so here "sym"
                ## In the recursive calls, HnH_info will be at 1, but the original value is kept for the rest of this section
                ## The same is true for symmetry_info
                symmetry_info_local = -sgn

                # feasibility check: without all stem vectors, one cannot be sure the system is R^{n+1}-feasible
                @suppress begin
                    tick()
                end
                d, lambda = isf_feas_HnH!(svec[1:nv]' .* Vt[:,T], -svTd*Vt[:,ivp], info, options, values)
                info.cput_lp += tok()
                info.flag > 0 && return 

                if any(d .!= Inf)       # system is feasible

                    if p > nvp
                        isf_rec_HnH!(Vt, svec, perm, d, 1, symmetry_info_local, info, options, values)
                        info.flag > 0 && return 
                    else
                        # no need to test symmetry because it's given right above
                        info.ns += 1
                        options.s && isf_storing_HnH!(svec, p, perm, symmetry_info_local, info)
                    end
                else

                    if options.sv == 2      # dual variable becomes a stem vector
                        # already in "sym"-covered status: if a new vector is useful it is necessary for Vt
                        isf_sv_add_HnH!(Vt, svec .* lambda[1:nvp], "asym", perm, info, options)
                    end
                end
            end

        else        # a LO problem is solved to check [svec;-1]
            @suppress begin
                tick()
            end
            
            d, lambda = isf_feas_HnH!(svec[1:nv]' .* Vt[1:n,T], -svTd*Vt[1:n,ivp], info, options, values)
            info.cput_lp += tok()
            info.flag > 0 && return 
            if any(d .!= Inf)      # system is feasible

                if p > nvp
                    isf_rec_HnH!(Vt, svec, perm, d, HnH_info, symmetry_info, info, options, values)
                    info.flag > 0 && return 
                else
                    info.ns += 2-HnH_info
                    options.s && isf_storing_HnH!(svec, p, perm, symmetry_info, info)
                end
            else

                if options.sv == 2      # dual variable becomes a stem vector

                    if HnH_info == 1
                        isf_sv_add_HnH!(Vt, svec.* lambda[1:nvp], "asym", perm, info, options)
                    else
                        if abs((svec.* lambda[1:nvp])' * Vt[n+1,[T;ivp]]) < options.tol_nonzero_q
                            isf_sv_add_HnH!(Vt, svec.* lambda[1:nvp], "asym", perm, info, options) ### the opposite thing compared to the nH version
                        else
                            isf_sv_add_HnH!(Vt, svec.* lambda[1:nvp], "sym", perm, info, options)
                        end
                    end
                    # here as the recursion can continue, one needs to modify manually the svsprod variables
                    if options.svsprod      # update of the recursive computation
                        if info.nb_stems_sym > length(svsprod_common_sym)   # stem added, so current size should be corrected
                            @suppress begin
                                tick()
                            end
                            tempo = zeros(Int,p)
                            tempo[perm[1:nv]] = svec[1:nv]
                            info.svsprod_sym = [info.svsprod_sym ; info.stems_sym[:,length(svsprod_common_sym)+1:info.nb_stems_sym]'* tempo]
                            # svsprod_common_sym = [svsprod_common_sym ; info.stems_sym[:,length(svsprod_common_sym)+1:info.nb_stems_sym]'* tempo]  # here, THIS should be removed
                            info.cput_cover += tok()
                        end
                        if info.nb_stems_asym > length(svsprod_common_asym)   # stem added, so current size should be corrected
                            @suppress begin
                                tick()
                            end
                            tempo = zeros(Int,p)
                            tempo[perm[1:nv]] = svec[1:nv]
                            info.svsprod_asym = [info.svsprod_asym ; info.stems_asym[:,length(svsprod_common_asym)+1:info.nb_stems_asym]'* tempo]
                            # svsprod_common_asym = [svsprod_common_asym ; info.stems_asym[:,length(svsprod_common_asym)+1:info.nb_stems_asym]' * tempo]    # here, THIS should be removed
                            info.cput_cover += tok()
                        end
                    end

                end
                # if HnH_info == 0 and infeasible, it might still be feasible in R^{n+1} (thus HnH_info -> 1)
                if HnH_info == 0 # [:,...] because we know it's the upper dimension
                    @suppress begin
                        tick()
                    end
                    d, lambda = isf_feas_HnH!(svec[1:nv]' .* Vt[:,T], -sign(vTd)*Vt[:,ivp], info, options, values)
                    info.cput_lp += tok()
                    info.flag > 0 && return 
                    if any(d .!= Inf)      # system is feasible

                        if p > nvp
                            isf_rec_HnH!(Vt, svec, perm, d, 1, -Int(sign(d[n+1])), info, options, values)
                            info.flag > 0 && return 
                        else
                            info.ns += 1
                            options.s && isf_storing_HnH!(svec, p, perm, -Int(sign(d[n+1])), info)
                        end
                    else

                        if options.sv == 2      # dual variable becomes a stem vector
                            isf_sv_add_HnH!(Vt, svec .* lambda[1:nvp], "asym", perm, info, options)
                        end
                    end
                end
            end
        end
    end     # of the vTd != 0 test
end