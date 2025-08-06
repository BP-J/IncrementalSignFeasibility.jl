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
`isf_first_HnH!(V, info, options, values)`
Modifies info and prepare the recursive call in the mixed case. 
Does the rank computation, the eventual stem vectors computations. 
Also contains the recursive calls and binary manipulations. 
"""
function isf_first_HnH!(Vt::Matrix, info::Info, options::Options, values::Values)

    p0 = size(Vt,2)
    if options.rational
        Vt, colsel, colout = isf_noncolin_r!(Vt)
    else
        Vt, colsel, colout = isf_noncolin!(Vt, options)
    end
    # now the vectors & right hand sides are normed
    
    p = size(Vt,2)
    n = size(Vt,1) - 1

    # Rank computation: if isf_first is used, options.algorithm should not be 0, 
    # so all algorithms use the QR factorization and no test required
    # Julia's QR factorization returns a complicated object with fields:
    # .Q, .R, .p, .P; p is the permutation vector and P the permutation matrix
    # if qr(A), one has A * F.P = F.Q * F.R

    F = qr(Vt[1:n,:], ColumnNorm()) # the second argument is the default (?) argument so there is a permutation
    Q = F.Q
    R = F.R
    perm = F.p
    
    if options.rational
        r = rankx(Vt[1:n,:])
        info.r = r
    else
        r = size(R, 1)
        while r > 0 && (norm(R[r,:]) < options.tol_nonzero_q)
            r = r-1
        end
        info.r = r
    end

    ### computation of all the stem vectors
    info.stems_sym = Matrix{Int64}(undef, p, 0)
    info.stems_asym = Matrix{Int64}(undef, p, 0)
    info.stems_sym_init = ElasticMatrix{Int64}(undef, p, 0)
    info.stems_asym_init = ElasticMatrix{Int64}(undef, p, 0)
    
    if options.sv >= 3

        isf_choose_wechelon_option!(Vt, Q, R, info, options)
        @suppress begin
            tick()
        end
        if options.wechelon
            if options.rational
                isf_allsv_HnH_Wechelon_r!(Vt, perm, info, options, values)
            else
                if r < n            # here the rank != the dimension, so it is better to look at the matrix R to avoid numerical issues
                    isf_allsv_HnH_Wechelon!([R[1:r, invperm(perm)]; Vt[n+1,:]], perm, info, options, values)
                else
                    isf_allsv_HnH_Wechelon!(Vt, perm, info, options, values)
                end
            end
        else
            if options.rational
                isf_allsv_HnH_r!(Vt, perm, info, options, values)
            else
                if r < n            # here the rank != the dimension, so it is better to look at the matrix R to avoid numerical issues
                    isf_allsv_HnH!([R[1:r, invperm(perm)]; Vt[n+1,:]], perm, info, options, values)
                else
                    isf_allsv_HnH!(Vt, perm, info, options, values)
                end
            end
        end
        info.cput_sv += tok()

    elseif options.sv > 0
        ### This is when some stem vectors are computed from the QR factorization
        if p > r # p = r means there is nothing to do
            @suppress begin
                tick()
            end
            if r < n
                isf_somesv_HnH!([R[1:r, invperm(perm)]; Vt[n+1,:]], perm, info, options)
            else
                isf_somesv_HnH!(Vt, perm, info, options)
            end
            info.cput_sv += tok()
        end
    end

    if mod(options.algorithm, 8) >= 6
        isf_choose_svsprods_option!(Vt, info, options)
    end

    # no need to get the 'common point' in this version

    svec = ones(Int,r)          # same convention as the H case
    Rinv = inv(R[1:r, 1:r])'    # used to compute the starting points

    if !options.withd
        d = []
    end

    while true
        # Associated direction
        if options.withd 
            d = Q[:,1:r] * (Rinv * svec)
            d = d / norm(d)
        end

        # recursive call
        # HnH_info initialized at 0: the recurrence of each of the 2^r starts in the normal state
        # symmetry_info as well    : initially the right part belonging to the asymmetric set is unknown

        if r < p 
            if options.svsprod  # prepare the recursive computation of the stem vector * current sign vector product
                @suppress begin
                    tick()
                end
                info.svsprod_sym  = info.stems_sym[perm[1:r],:]' * svec
                info.svsprod_asym = info.stems_asym[perm[1:r],:]' * svec
                info.cput_cover += tok()
            end 
            if options.withd 
                isf_rec_HnH!(Vt, svec, perm, d, 0, 0, info, options, values)
            else                # algorithm 15 := AD4, only stem vectors
                isf_rec_nod_HnH!(Vt, svec, perm, 0, 0, info, options, values)    
            end
            info.flag > 0 && return info
        else
            info.ns += 1
            isf_storing!(svec, p, perm, info)
        end

        # next binary vector

        s = isf_decrement_s(svec[2:r])
        if s == []
            break   # for the while true
        end
        svec[2:r] = s

    end

#-------------------------------------------------------------------------------------------------------------------------------
# if some columns were colinear, they have been removed so they must be added again
#-------------------------------------------------------------------------------------------------------------------------------

    if p0 != p
        eliminated_columns = setdiff(1:p0, colsel)
        mat_sign_vec = zeros(Int, info.ns, p0)
        mat_sign_vec[:, colsel] = hcat(info.s...)'  # puts a list of vectors into a matrix (same 'dimensions')
        # now the 'matrix' of signs vectors must be completed for the removed indices, what was computed previously is put at the correct lines
        # then for each index that was removed, either the coordinate is duplicated or duplicated then * -1
        for elim in eliminated_columns
            if colout[elim] > 0
                mat_sign_vec[:,elim] =  mat_sign_vec[:,abs(colout[elim])]
            else
                mat_sign_vec[:,elim] = -mat_sign_vec[:,abs(colout[elim])]
            end
        end
        info.s = [mat_sign_vec[i,:] for i in 1:info.ns]

        ## treating the stem vectors
        if options.sv > 0

            # preliminary work - to put at the correct place the number of times a vector is duplicated
            info_colinearity = ones(Int, p)
            for elim in colout[eliminated_columns]
                info_colinearity[findfirst(t -> t == abs(elim), colsel)] += 1
            end
            # after that, for every i in [1:p], info_colinearity[i] gives the number of vectors colinear to v_i (including itself)

            STEMS = info.stems_sym
            info.stems_sym = zeros(Int, p0, info.nb_stems_sym)
            info.stems_sym[colsel, :] = STEMS
            # now the lines of removed indices are empty
            # but stem vectors must be added, by swapping the coordinates at indices of colinear vectors
            # it is done in the following way:
            # for each stem vector, for each index, it is verified if it is an index to duplicate
            #       for example, if coordinate i1 of the current stem is duplicated 3 times, and i2 2 times 
            #       [the coordinates for these indices being ±1], then the current stem must be duplicated 3*2 = 6 times
            #       which leads to 5 new vectors - 3 with coordinate i2 unchanged and the different positions for i1, 
            #       and 3 with coordinate i2 changed and the different positions for i1; the one with no changes = the current stem

            for N in 1:info.nb_stems_sym
                ajout_N = info.stems_sym[:,N]   # this one to duplicate, not itself added at the end (already there)
                base_size = 1 # number of times to be duplicated

                for index_1 in 1:p
                    # for every index_1, if the coordinate !=0 (otherwise no need to duplicate) and v_index_1 has a colinear one (otherwise no need to duplicate)
                    if info_colinearity[index_1] >= 2 && info.stems_sym[colsel[index_1],N] != 0  # if == 1, vector has no colinear ones
                        ajout_N = repeat(ajout_N, inner=(1,1), outer=(1,info_colinearity[index_1])) # duplication
                        # searching the indices that are duplicates of V[:,colsel[index_1]]
                        indices = eliminated_columns[findall(t -> abs(t) == colsel[index_1], colout[eliminated_columns])]
                        # changing the submatrices of ajout_N
                        for index_2 in 1:info_colinearity[index_1]-1# this number of changes to do
                            ajout_N[indices[index_2], base_size*(index_2)+1:base_size*(index_2+1)] = sign(colout[indices[index_2]]) * ajout_N[colsel[index_1], base_size*(index_2)+1:base_size*(index_2+1)]
                            ajout_N[colsel[index_1], base_size*(index_2)+1:base_size*(index_2+1)] *= 0
                        end
                        base_size *= info_colinearity[index_1]
                    end
                end
                info.stems_sym = [info.stems_sym ajout_N[:,2:base_size]]
            end
            info.nb_stems_sym = size(info.stems_sym, 2)

            # for the asymmetric ones

            STEMS = info.stems_asym
            info.stems_asym = zeros(Int, p0, info.nb_stems_asym)
            info.stems_asym[colsel, :] = STEMS
            number_initial = info.nb_stems_asym

            # the bloc for the pairwise ones (by convention, they are "asym" here)
            # as some vectors are colinear, this means alpha_i v_i = alpha_j v_j, so there is a stem with coordinates just between indices i and j
            # this can be easily checked and added. What is trickier is, if for instance vectors 3,6,7 are colinear, that there is a stem between 3 & 6, 3 & 7 but also 6 & 7
            # therefore there is a tedious verification to add the correct (number of) stem vectors

            Temp = zeros(Int, p0, 0)
            for index_1 in 1:p                              # reference to info_colinearity
                if info_colinearity[index_1] > 1            # otherwise nothing to do
                    indices = findall(t -> t == index_1, abs.(colout))
                    N = length(indices)                     # >= 2 due to info_colinearity[index_1] > 1
                    for index_2 in 1:N-1                    # double loop on
                        for index_3 in index_2+1:N          # the indices
                            STEM = zeros(Int, p0)           # creation of the new stem
                            STEM[indices[index_2]] = +sign(colout[indices[index_2]])
                            STEM[indices[index_3]] = -sign(colout[indices[index_3]])
                            info.nb_stems_asym += 1
                            Temp = [Temp STEM]
                        end
                    end
                end
            end
            info.stems_asym = [info.stems_asym Temp]

            for N in 1:number_initial
                ajout_N = info.stems_asym[:,N]   # this one to duplicate, not itself added at the end (already there)
                base_size = 1 # number of times to be duplicated

                for index_1 in 1:p
                    if info_colinearity[index_1] >= 2 && info.stems_asym[colsel[index_1],N] != 0  # if == 1, vector has no colinear ones
                        ajout_N = repeat(ajout_N, inner=(1,1), outer=(1,info_colinearity[index_1])) # duplication
                        # searching the indices that are duplicates of V[:,colsel[index_1]]
                        indices = eliminated_columns[findall(t -> abs(t) == colsel[index_1], colout[eliminated_columns])]
                        # changing the submatrices of ajout_N
                        for index_2 in 1:info_colinearity[index_1]-1# this number of changes to do
                            ajout_N[indices[index_2], base_size*(index_2)+1:base_size*(index_2+1)] = sign(colout[indices[index_2]]) * ajout_N[colsel[index_1], base_size*(index_2)+1:base_size*(index_2+1)]
                            ajout_N[colsel[index_1], base_size*(index_2)+1:base_size*(index_2+1)] *= 0
                        end
                        base_size *= info_colinearity[index_1]
                    end
                end
                info.stems_asym = [info.stems_asym ajout_N[:,2:base_size]]
            end
            
            info.nb_stems_asym = size(info.stems_asym,2)

            info.stem_sizes_sym  = vec(sum(abs.(info.stems_sym), dims=1))                               # the sum returns a 1 x [size] matrix, not a vector
            info.stem_sizes_asym = vec(sum(abs.(info.stems_asym), dims=1))                              # the sum returns a 1 x [size] matrix, not a vector
        end
    end

end