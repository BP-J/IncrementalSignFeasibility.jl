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
`isf_rec_nod!(V, svec, perm, type, info, options, values)`

Recursive procedures forming the S-tree, without using directions. 
These versions with _nod are shorter and 'simpler' than the versions 
with, as they 'only' look if a current sign vector covers one of the
stem vectors. 
The H and nH versions are extremely similar, the main differences being:
the 2 asymmetric starts (but not apparent here) and the asymmetry of the 
test, which is done the `isf_cover_stem_*` functions. 

The mixed HnH version is different enough that it was done in another
function, though it remains quite similar. 
"""
function isf_rec_nod!(V::Matrix, svec::Vector{Int}, perm::Vector{Int}, type::String, info::Info, options::Options, values::Values)

    info.flag = values.success

    p   = size(V,2)
    r   = info.r
    nv  = length(svec)
    nvp = nv+1

    ## no heuristic about vector order in these versions

#-------------------------------------------------------------------------------------------------------------------------------
# Store the 'current' svsprod that will be modified for the recursive calls
#-------------------------------------------------------------------------------------------------------------------------------

    if options.svsprod              # it is necessary to keep the 'original' value in the node; in the first subbranch the value changes
        @suppress begin
            tick() 
        end
        svsprod_common_sym  = info.svsprod_sym
        svsprod_common_asym = info.svsprod_asym
        info.cput_cover += tok()
    end

    # first sign tested
    svec = [svec; -1]
    next_str = "-1"

    info.nb_sv_checks += 1
    if type == "H"
        if options.svsprod
            @suppress begin
                tick()
            end
            info.svsprod_sym  -= info.stems_sym[perm[nvp],:]
            isc = isf_cover_stem_H_bis(svec, info, options)                         # covering test, recursive version
            info.cput_cover += tok()
        else
            @suppress begin
                tick()
            end
            isc = isf_cover_stem_H(svec, perm, info, options)                       # covering test, full version
            info.cput_cover += tok()
        end
    else
        if options.svsprod
            @suppress begin
                tick()
            end
            info.svsprod_sym  -= info.stems_sym[perm[nvp],:]
            info.svsprod_asym -= info.stems_asym[perm[nvp],:]
            isc = isf_cover_stem_nH_bis(svec, info, options)                        # covering test, recursive version
            info.cput_cover += tok()
        else
            @suppress begin
                tick()
            end
            isc = isf_cover_stem_nH(svec, perm, info, options)                      # covering test, full version
            info.cput_cover += tok()
        end
    end

    if isc#!isempty(isc)        # isc != empty, so there is an index covered: current svec is infeasible
        info.nb_sv_detect += 1
        ## recursion stopped for this one...

        ## ... but then the other sign is necessarily feasible
        svec[nvp] = +1

        # recursive call

        if p > nvp
            if type == "H"
                if options.svsprod
                    @suppress begin
                        tick()
                    end
                    info.svsprod_sym  = svsprod_common_sym + vec(info.stems_sym[perm[nvp],:])
                    info.cput_cover += tok()
                end
            else    # "nH"
                if options.svsprod
                    @suppress begin
                        tick()
                    end
                    info.svsprod_sym  = svsprod_common_sym + info.stems_sym[perm[nvp],:]
                    info.svsprod_asym = svsprod_common_asym + info.stems_asym[perm[nvp],:]
                    info.cput_cover += tok()
                end
            end
            isf_rec_nod!(V, svec, perm, type, info, options, values)

            info.flag > 0 && return info
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end
    else
        ## current one is feasible, so recursion continued; the other one is tested after

        # recursive call

        if p > nvp
            isf_rec_nod!(V, svec, perm, type, info, options, values)
            info.flag > 0 && return info
        else
            info.ns += 1
            options.s && isf_storing!(svec, p, perm, info)
        end

        svec[nvp] = +1
        next_str = "+1"

        info.nb_sv_checks += 1
        if type == "H"
            if options.svsprod
                @suppress begin
                    tick()
                end
                info.svsprod_sym = svsprod_common_sym + info.stems_sym[perm[nvp],:]
                isc = isf_cover_stem_H_bis(svec, info, options)
                info.cput_cover += tok()
            else
                @suppress begin
                    tick()
                end
                isc = isf_cover_stem_H(svec, perm, info, options)
                info.cput_cover += tok()
            end
        else
            if options.svsprod
                @suppress begin
                    tick()
                end
                info.svsprod_sym  = svsprod_common_sym + info.stems_sym[perm[nvp],:]
                info.svsprod_asym = svsprod_common_asym + info.stems_asym[perm[nvp],:]
                isc = isf_cover_stem_nH_bis(svec, info, options)
                info.cput_cover += tok()
            else
                @suppress begin
                    tick()
                end
                isc = isf_cover_stem_nH(svec, perm, info, options)
                info.cput_cover += tok()
            end
        end

        if !isc#isempty(isc) # other one also feasible
            info.nb_sv_detect += 1

            if p > nvp
                isf_rec_nod!(V, svec, perm, type, info, options, values)
                
                info.flag > 0 && return info
            else
                info.ns += 1
                options.s && isf_storing!(svec, p, perm, info)
            end
        end
    end
end

#-----------------------------------------------------------------------

"""
`isf_rec_nod_HnH!(V, svec, perm, type, info, options, values)`

Recursive procedures forming the S-tree, without using directions. 
These versions with _nod are shorter and 'simpler' than the versions 
with, as they 'only' look if a current sign vector covers one of the
stem vectors. 
The H and nH versions are extremely similar, but a  bit simpler so 
done together in a separate function. 

The mixed HnH version is different enough that it was done in another
function, though it remains quite similar. It uses additional inputs
HnH_info, symmetry_info, which tell the current state of the recursion. 
- HnH_info == 0 means the recursion is in the 'natural' state, (R^n)
- HnH_info == 1 means the recursion is in the 'upgraded' state, i.e.,
  the current branch is V-infeasible but Vt-feasible
symmetry_info is a sign telling which sign vector is really in S
(one could use only symmetry_info and "if symmetry_info != 0"...)
"""
function isf_rec_nod_HnH!(Vt::Matrix, svec::Vector{Int}, perm::Vector{Int}, HnH_info::Int, symmetry_info::Int, info::Info, options::Options, values::Values)
    info.flag = values.success

    p   = size(Vt,2)
    nv  = length(svec)
    nvp = nv+1

#-------------------------------------------------------------------------------------------------------------------------------
# Store the 'current' svsprod that will be modified for the recursive calls
#-------------------------------------------------------------------------------------------------------------------------------

    @suppress begin
        tick() 
    end
    if options.svsprod
        svsprod_common_sym  = info.svsprod_sym
        svsprod_common_asym = info.svsprod_asym
    end
    info.cput_cover += tok()

    ## no heuristic about vector order in these versions

    svec_base = copy(svec)
    svec = [svec_base;+1]
    next_str = "+1"
    info.nb_sv_checks += 1

    @suppress begin
        tick()
    end
    if options.svsprod
        info.svsprod_asym += info.stems_asym[perm[nvp],:]

        # modifying the symmetric one only if basic state
        if HnH_info == 0
            info.svsprod_sym += info.stems_sym[perm[nvp],:]
        end
        covering_index, detail, sgn = isf_cover_stem_HnH_bis(svec, HnH_info, info, options)
    else
        covering_index, detail, sgn = isf_cover_stem_HnH(svec, perm, HnH_info, info, options)
    end
    info.cput_cover += tok()

    if covering_index # covering_index is not [], so is the index of a covered stem vector (there might be much more)
        ## this means the current svec is infeasible... (*)
        info.nb_sv_detect += 1
        if detail == "asym"     # means the current sign vector is infeasible even in dimension n+1, so the recursion is stopped

            # recursion is fully stopped, nothing to do

        else        # covering_index != [] means detail ∈ ["sym", "asym"], so here "sym"

            symmetry_info_local = -sgn
            if p > nvp 
                isf_rec_nod_HnH!(Vt, svec, perm, 1, symmetry_info_local, info, options, values)
                info.flag > 0 && return info
            else
                info.ns += 1        # here we know the current level has been stopped AND is the final level, so just that is fine
                if options.s
                    push!(info.s, zeros(Int, p))
                    if symmetry_info_local == -1
                        info.s[info.ns][perm] = - svec # update to put in the right one
                    else
                        info.s[info.ns][perm] = svec # update to put in the right one
                    end
                end
            end

        end # of detail == ...

        ## (*) but if the first tested sign (0) is infeasible, then the other is necessarily feasible

        svec = [svec_base;-1] # the other sign
        next_str = "-1"

        if p > nvp 
            if options.svsprod
                info.svsprod_asym = svsprod_common_asym - info.stems_asym[perm[nvp],:]

                if HnH_info == 0
                    info.svsprod_sym  = svsprod_common_sym - info.stems_sym[perm[nvp],:]
                end
            end
            isf_rec_nod_HnH!(Vt, svec, perm, HnH_info, symmetry_info, info, options, values)
            info.flag > 0 && return info
        else
            if symmetry_info == 0   # the recurrence is still symmetric so both are added
                info.ns += 2
                if options.s
                    push!(info.s, zeros(Int,p))
                    info.s[info.ns-1][perm] = svec 
                    push!(info.s, zeros(Int,p))
                    info.s[info.ns][perm] = - svec 
                end
            else # symmetry_info ∈ {±1}; the recurrence is not symmetric anymore so only the right one is added
                info.ns += 1
                if options.s
                    push!(info.s, zeros(Int,p))
                    if symmetry_info == -1
                        info.s[info.ns][perm] = - svec # update to put in the right one
                    else
                        info.s[info.ns][perm] = svec # update to put in the right one
                    end
                end
            end
        end
    else    # !isempty(covering_index) is false, so covering_index is empty: the recurrence is continued without changing the status
        if p > nvp
            isf_rec_nod_HnH!(Vt, svec, perm, HnH_info, symmetry_info, info, options, values)
            info.flag > 0 && return info
        else
            if symmetry_info == 0   # the recurrence is still symmetric so both are added
                info.ns += 2
                if options.s
                    push!(info.s, zeros(Int,p))
                    info.s[info.ns-1][perm] = svec 
                    push!(info.s, zeros(Int,p))
                    info.s[info.ns][perm] = - svec 
                end
            else # symmetry_info ∈ {±1}; the recurrence is not symmetric anymore so only the right one is added
                info.ns += 1
                if options.s
                    push!(info.s, zeros(Int,p))
                    if symmetry_info == -1
                        info.s[info.ns][perm] = - svec # update to put in the right one
                    else
                        info.s[info.ns][perm] = svec # update to put in the right one
                    end
                end
            end
        end

        # now the other sign is checked
        svec = [svec_base;-1]
        info.nb_sv_checks += 1
        next_str = "-1"

        @suppress begin
            tick()
        end
        if options.svsprod
            info.svsprod_asym = svsprod_common_asym - info.stems_asym[perm[nvp],:]
            if HnH_info == 0
                info.svsprod_sym = svsprod_common_sym - info.stems_sym[perm[nvp],:]
            end
            covering_index, detail, sgn = isf_cover_stem_HnH_bis(svec, HnH_info, info, options)
        else
            covering_index, detail, sgn = isf_cover_stem_HnH(svec, perm, HnH_info, info, options)
        end
        info.cput_cover += tok()

        if !covering_index # the recurrence is continued
            info.nb_sv_detect += 1
            if p > nvp
                isf_rec_nod_HnH!(Vt, svec, perm, HnH_info, symmetry_info, info, options, values)
                info.flag > 0 && return info
            else
                if symmetry_info == 0   # the recurrence is still symmetric so both are added
                    info.ns += 2
                    if options.s
                        push!(info.s, zeros(Int,p))
                        info.s[info.ns-1][perm] = svec 
                        push!(info.s, zeros(Int,p))
                        info.s[info.ns][perm] = - svec 
                    end
                else # symmetry_info ∈ {±1}; the recurrence is not symmetric anymore so only the right one is added
                    info.ns += 1
                    if options.s
                        push!(info.s, zeros(Int,p))
                        if symmetry_info == -1
                            info.s[info.ns][perm] = - svec # update to put in the right one
                        else
                            info.s[info.ns][perm] = svec # update to put in the right one
                        end
                    end
                end
            end

        else
            if detail == "asym" # means the current sign vector is infeasible even in dimension n+1, so the recursion is stopped

                # recursion is fully stopped, nothing to do

            else        # covering_index != [] means detail ∈ ["sym", "asym"], so here "sym"

                symmetry_info_local = -sgn
                if p > nvp 
                    isf_rec_nod_HnH!(Vt, svec, perm, 1, symmetry_info_local, info, options, values)
                    info.flag > 0 && return info
                else
                    info.ns += 1        # here we know the current level has been stopped AND is the final level, so just that is fine
                    if options.s
                        push!(info.s, zeros(Int, p))
                        if symmetry_info_local == -1
                            info.s[info.ns][perm] = - svec # update to put in the right one
                        else
                            info.s[info.ns][perm] = svec # update to put in the right one
                        end
                    end
                end
    
            end # of detail == ...
        end ## of the second isempty(covering_index)
    end ## of the first !isempty(covering_index)
end