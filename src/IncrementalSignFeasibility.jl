module IncrementalSignFeasibility

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

# useful packages
using Dates                             # to print the date
using LinearAlgebra                     # null space, qr factorization, norms...
using LinearAlgebraX                    # exact computations for the raitonal case
using Gurobi                            # Gurobi solver
using GLPK                              # GLPK solver
using JuMP                              # Mathematical Programming package
using Printf                            # printing with @printf
using TickTock                          # simple time measurement
using Suppressor                        # removing undesirable printings 
using ElasticArrays                     # for a new structure
using RowEchelon                        # for the alternate stem computation with echelon form
using StatsBase                         # for the sample function (generation of data)
using Random                            # for randomness in the generation
using Combinatorics                     # for some instance generation
using DataFrames                        # to gather data in some organized output
using DelimitedFiles                    # to treat various types of files

### other files
include("isf_tools.jl")                     # short functions/structures/... used in the rest
# has to be done first as it defines types and so on
include("isf_allsv_from_indices_H.jl")      # individual computation for a subset of columns
include("isf_allsv_from_indices_nH.jl")     # individual computation for a subset of columns
include("isf_allsv_from_indices_HnH.jl")    # individual computation for a subset of columns
include("isf_allsv_H.jl")                   # stem vectors in the homogeneous case (V)
include("isf_allsv_H_Wechelon.jl")          # stem vectors in the homogeneous case (V) with echelon forms 
include("isf_allsv_H_Wechelon_rec.jl")      # stem vectors in the homogeneous case (V) with echelon forms - rec part
include("isf_allsv_nH.jl")                  # stem vectors in the non-homogeneous case (Vt)
include("isf_allsv_nH_Wechelon.jl")         # stem vectors in the non-homogeneous case (Vt) with echelon forms
include("isf_allsv_nH_Wechelon_rec.jl")     # stem vectors in the non-homogeneous case (Vt) with echelon forms - rec part
include("isf_allsv_HnH.jl")                 # stem vectors in the hybrid case (V - Vt)
include("isf_allsv_HnH_Wechelon.jl")        # stem vectors in the hybrid case (V - Vt) with echelon forms - rec part
include("isf_allsv_HnH_Wechelon_rec.jl")    # stem vectors in the hybrid case (V - Vt) with echelon forms - rec part
#include("isf_benchmark_values.jl")          # tool used in other benchmarking functions
include("isf_bin_operations.jl")            # manipulations on binary vectors
include("isf_binary_tree_update.jl")        # manipulations on binary vectors
include("isf_bounds.jl")                    # Winder / Zaslavsky bounds
include("isf_check.jl")                     # to check options at the beginning
include("isf_choose_hard_options.jl")       # modifies the difficult options
include("isf_coherence.jl")                 # to check coherence of options
include("isf_cover_stem.jl")                # covering test
include("isf_cover_stem_bis.jl")            # covering test - second version
include("isf_echelonning.jl")               # echelonned form 
include("isf_feas.jl")                      # feasibility problems
include("isf_first_H.jl")                   # initialization
include("isf_first_nH.jl")                  # initialization
include("isf_first_HnH.jl")                 # initialization
include("isf_first_rc2018.jl")              # initialization
include("isf_generation_affine.jl")         # creating affine instances
include("isf_generation_linear.jl")         # creating linear instances
include("isf_nod_lazy_completion.jl")       # recursive completion of the stems
include("isf_nod_lazy.jl")                  # lazy AD5 method
include("isf_noncolin.jl")                  # removing the duplicate vectors in the data
include("isf_ntype.jl")                     # heuristic B
include("isf_rec_H.jl")                     # recursive call
include("isf_rec_nH.jl")                    # recursive call
include("isf_rec_HnH.jl")                   # recursive call
include("isf_rec_nod.jl")                   # recursive call (stem vectors)
include("isf_rec_rc2018.jl")                # recursive call 
include("isf_somesv_H.jl")                  # a few stem vectors
include("isf_somesv_nH.jl")                 # a few stem vectors
include("isf_somesv_HnH.jl")                # a few stem vectors
include("isf_storing.jl")                   # when length is p, recursion stops
include("isf_sv_add.jl")                    # converts dual variables into stem vectors

### exports
export isf, readdlm, Options, Info, Values, isf_get_options, isf_get_info, isf_get_values, options_from_algo, options_rationals # not the end of isf_tools.jl
export isf_allsv_H!, isf_allsv_from_indices_H!, isf_allsv_from_indices_H_r!, isf_allsv_nH!, isf_allsv_from_indices_nH!, isf_allsv_from_indices_nH_r!, isf_allsv_HnH!, isf_allsv_from_indices_HnH!, isf_allsv_from_indices_HnH_r!
export isf_allsv_H_Wechelon!, isf_allsv_H_Wechelon_rec!, isf_allsv_H_Wechelon_rec_r!, isf_allsv_nH_Wechelon!, isf_allsv_nH_Wechelon_rec!, isf_allsv_nH_Wechelon_rec_r!, isf_allsv_HnH_Wechelon!, isf_allsv_HnH_Wechelon_rec!, isf_allsv_HnH_Wechelon_rec_r!
export isf_binary_minus, isf_binary_plus, isf_decrement_s, isf_increment_s, isf_binary_tree_update, isf_central_max, isf_noncentral_max, isf_check, isf_choose_wechelon_option!, isf_choose_svsprods_option!, isf_coherence!
export isf_cover_stem_H, isf_cover_stem_nH, isf_cover_stem_HnH, isf_cover_stem_H_bis, isf_cover_stem_nH_bis, isf_cover_stem_HnH_bis, echelonned_form!, isf_feas!, isf_feas_HnH!, isf_noncolin!, isf_noncolin_r!, isf_ntype_H, isf_ntype_nH, isf_ntype_HnH, isf_nod_lazy_H!, isf_nod_lazy_nH!, isf_nod_completion_H!, isf_nod_completion_nH_sym!, isf_nod_completion_nH_asym!
export isf_first_H!, isf_first_nH!, isf_first_HnH!, isf_first_rc2018!, isf_rec_rc2018!, isf_rec_H!, isf_rec_nH!, isf_rec_HnH!, isf_rec_nod!, isf_rec_nod_HnH!
export isf_somesv_H!, isf_somesv_nH!, isf_somesv_HnH!, isf_storing!, isf_storing_HnH!, isf_sv_add_H!, isf_sv_add_nH!, isf_sv_add_HnH!

# --------------------------------------------------

###

# This file sets up the useful tools used for the computation of the 
# 'sign vector set' of a matrix Vt = [V;T] representing an arrangement 
# of hyperplanes. Vt is assumed to be of the form [vt_1, ... vt_p]
# where vt_i = (v_i, τ_i) is of dimension n+1, the last 'coordinate' playing 
# particuar role, the hyperplane H_i being defined by {x ; (v_i, x) = τ_i}
#
# The goal is to compute the sign vectors s in {-1, +1}^p such that 
# the system s .* (V' * x - T') > 0 has a solution xs of size n. 
#
# Important notation: it is possible to compute a central arrangement,
# i.e., where all the τ_i = 0 (so T = 0).The code checks whether 
# this holds and if so launches slightly simpler computations. 

# Various options allow the user to tune the code, though it is assumed 
# the options are coherent with the data: for instance, if the arrangement
# is not symmetric, it is assumed the user did not say it was.

# TBD

# . options.symmetry: [boolean] describes if the arrangement is symmetric: 
# if unknown, say false and the code checks it
# . options.s: [integer] if 1, store the feasible s's in info.s
# . options.bestv (integer in [0:3]): can be used to modify the order in
# which the vectors (columns of V) are considered for constructing the
# S-tree (only the value 3 is actually documented and evaluated in the
# paper):
# == 0: no modification of the order;
# == 3: a reodering of the vectors is made so that the number of nodes
#       of the S-tree decreases (and therefore the number of LOP to
#       solve);
# . options.dvnear0 (logical): is true, the S-tree algorithm constructs
#   two descendants when v'*d is in a specific interval surrounding zero
#   (v is the new considered vector)

# . options.agorithm: [integer] as described, except for D4, 
# other letters suppose the modifications with previous letters 
# are also included: B means A and B, D2 means A, B, C, D2
# == 0:  original algorithm from Rada and Cerny [no modification]
# == 1:  RC algorithm with modification A  (At)
# == 2:  RC algorithm with modification B  (Bt)
# == 3:  RC algorithm with modification C  (Ct)
# == 4:  RC algorithm with modification D1 (Dt_1)
# == 5:  RC algorithm with modification D2 (Dt_2)
# == 6:  RC algorithm with modification D3 (Dt_3)
# == 7:  RC algorithm with modification D4 (Dt_4)
# == 8:  H-nH algorithm with no modification 
# == 9:  H-nH algorithm with modification A  (At)
# == 10: H-nH algorithm with modification B  (Bt)
# == 11: H-nH algorithm with modification C  (Ct)
# == 12: H-nH algorithm with modification D1 (Dt_1)
# == 13: H-nH algorithm with modification D2 (Dt_2)
# == 14: H-nH algorithm with modification D3 (Dt_3)
# == 15: H-nH algorithm with modification D4 (Dt_4)
#
# On return
# . info is a structure giving information on the run
# - info.flag: diagonis
# == 0: the required job has been realized;
# == 1: an input argument is wrong;
# == 2: nothing to do since V has no column,
# == 3: one of the vectors V(:,i) vanishes (then no solution),
# == 4: a pointed cone and its handle are not compatible
# (implementation or rounding error);
# == 5: at least a vector is opposite to another one, in which case
# the problem has no solution;
# == 6: the QP solver fails;
# == 9: a technical problem is encountered, which requires
# improvements of the code;
# - info.ns = half the number of sign vectors;
# - info.nsc = number of sign vectors in info.sclst;
# - info.s = half of the feasible sign vectors if options.s == true;
# - info.sclst = cell array of size (info.nsc,1) containing half of a
# root of infeasible sign vectors if options.sclst == 1; half of
# the infeasible sign vectors if options.sclst == 2; the rows of
# info.sclst is itself a cellarray of size (1,p) that gives binary
# vectors (when options.sclst == 1, some of the cell may by empty,
# meaning that they can take the values 0 or 1);

const GUROBI_ENV = Gurobi.Env()

"""
`info = isf(V, options)`
Realize the computations of the isf code. 
"""
function isf(Vt::Matrix, options::Options)
    ## main fuction, initial call
    ## calls the other (sub)functions along the warranty

    @suppress begin
        tick()
    end

    values = Values(0,1,2,3,4,5,6)    
    info = isf_get_info()

    if ~(options isa Options)
        info.flag = values.fail_argument
        throw(DomainError(options, "options is not well initialized")) 
        return info
    end

    options, info = isf_check(options, values, info)
    if info.flag > 1
        return info
    end

    ## slightly different depending on options.symmetry
    ## it is assumed that, if options.symmetry = true, the user
    ## checked it and doesn't use a matrix whose last row is the τ_i
    ## if options.symmetry = false, the last row is checked
    if !options.symmetry
        nn, p = size(Vt)
        n = nn-1
        Ttemp = Vt[nn, :] # the right-hand sides to be checked
        if norm(Ttemp) < 100*eps()
            options.symmetry == true
            Vt = Vt[1:n,:]          # the last row is nearly 0 and thus removed
        end
    else
        n, p = size(Vt)
    end

    # dimensions checked
    if p == 0
        info.ns = 0
        info.flag = values.fail_argument
        return info
    end

    if n == 0
        info.ns = 0
        info.flag = values.fail_argument
        return info
    end

    # incremental/recursive calls
    
    if options.algorithm == 0       # original RC, no modifications at all
        isf_first_rc2018!(Vt, info, options, values)
    elseif options.algorithm <= 7   # ISF with the modifications
        if options.symmetry
            isf_first_H!(Vt, info, options, values)
        else
            isf_first_nH!(Vt, info, options, values)
        end
    else                            # HnH algorithm, only if options.symmetry == false
        isf_first_HnH!(Vt, info, options, values)
    end

    info.cput_total = tok()

    if options.symmetry
        info.ns = 2*info.ns
    end

    return info
end

end
