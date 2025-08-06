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
# Not a particular function but defines some basic pieces reused later
# such as structures, type converting functions and so on.


#-----------------------------------------------------------------------
##### Structures
#-----------------------------------------------------------------------

mutable struct Options
    algorithm::Int
    bestv::Int
    dvnear0::Bool
    rational::Bool
    s::Bool
    # sclst::Bool
    sv::Int
    svsprod::Bool
    symmetry::Bool
    tol_coordinates::Float64
    tol_nonzero_q::Float64
    wechelon::Bool
    withd::Bool
end

mutable struct Values
    success::Int
    fail_argument::Int
    fail_zero_column::Int
    fail_solver::Int
    fail_solver_special::Int
    fail_not_positive::Int
    fail_technicality::Int
end

mutable struct Info
    flag::Int                                           # 0 if everything goes well, != 0 otherwise
    ns::Int                                             # number of sign vectors
    r::Int                                              # rank of V
    s::Vector{Vector{Int64}}                            # sign vectors
    # sclst::Set{SubArray{Bool, 1, Matrix{Bool}, Tuple{Base.Slice{Base.OneTo{Int64}}, Int64}, true}}    # sign vectors of the complementary set
    sclst::Set{Vector{Bool}}                            # sign vectors of the complementary set
    nb_losolve::Int64                                   # number of LO problems
    nb_feaslop::Int64                                   # number of   feasible LO problems
    nb_infeaslop::Int64                                 # number of infeasible LO problems
    nb_stems_sym::Int64                                 # number of  symmetric stem vectors
    nb_stems_asym::Int64                                # number of asymmetric stem vectors
    nb_sv_checks::Int64                                 # number of stem vectors tests
    nb_sv_detect::Int64                                 # number of stem vectors tests successful
    nb_duplicated_stems_sym::Int64                      # number of duplicated  symmetric stem vectors
    nb_duplicated_stems_asym::Int64                     # number of duplicated asymmetric stem vectors
    svsprod_sym::Vector{Int64}                          # recursive computation of the stem vector test,  symmetric part
    svsprod_asym::Vector{Int64}                         # recursive computation of the stem vector test, asymmetric part
    cput_lp::Float64                                    # time for the LO problems
    cput_sv::Float64                                    # time to obtain the stem vectors
    cput_cover::Float64                                 # time for the covering tests
    cput_total::Float64                                 # total time
    stems_sym::Matrix{Int64}                            # matrix of the  symmetric stem vectors
    stems_asym::Matrix{Int64}                           # matrix of the asymmetric stem vectors
    stems_sym_init::ElasticMatrix{Int64}                # precomputation of the stem vectors
    stems_asym_init::ElasticMatrix{Int64}               # precomputation of the stem vectors
    stem_sizes_sym::Vector{Int64}                       # sizes of the stem vectors
    stem_sizes_asym::Vector{Int64}                      # sizes of the stem vectors
    stem_zero_indices_sym::Vector{Vector{Int64}}        # certain ordering of the stem depending on the permutation
    stem_zero_indices_asym::Vector{Vector{Int64}}       # certain ordering of the stem depending on the permutation
end

#-----------------------------------------------------------------------
##### Basic initializations 
#-----------------------------------------------------------------------

"""
`options = isf_get_options()`
Initialization of options to default values. 
"""
function isf_get_options()
    ## creates the options from the matrix Vt alone
    options = Options(0, 0, false, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    return options
end

"""
`options = isf_get_options()`
Initialization of info to default values (modified by the algorithm). 
"""
function isf_get_info()
    info = Info(0, 0, 0, Vector{Vector{Int64}}(undef, 0), Set(), 0, 0, 0, 0, 0, 0, 0, 0, 0, Vector{Int64}(undef,0), Vector{Int64}(undef,0), 
    0.0, 0.0, 0.0, 0.0, Matrix{Int64}(undef, 0, 0), Matrix{Int64}(undef, 0, 0), ElasticMatrix{Int64}(undef, 0, 0), ElasticMatrix{Int64}(undef, 0, 0), Vector{Int64}(undef,0), Vector{Int64}(undef,0), 
    Vector{Vector{Int64}}(undef, 0), Vector{Vector{Int64}}(undef, 0))
    return info
end

"""
`options = isf_get_values()`
Initialization of values to default values (unused). 
"""
function isf_get_values()
    values = Values(0,1,2,3,4,5,6)
    return values
end

#-----------------------------------------------------------------------
##### Algorithmic initializations 
#-----------------------------------------------------------------------

"""
`options = options_from_algo(algo, symmetry, tol_coordinates, tol_nonzero)`
Adjustment of the options to the algorithm. 
"""
function options_from_algo(algo::Int, symmetry::Bool, tol_coordinates::Float64=100000*eps(), tol_nonzero::Float64=1000*eps())
    options = isf_get_options()
    options.symmetry = symmetry
    options.tol_coordinates = tol_coordinates
    options.tol_nonzero_q = tol_nonzero

    am8 = algo % 8
    if 2 <= am8 <= 6
        options.dvnear0 = true
        if 3 <= am8 <= 6
            options.bestv = 3
            options.sv = am8-3 # corresponds so just -3 instead of tests
        end
    end
    if am8 != 7
        options.withd = true
    else
        options.sv = 3
        options.withd = false
    end
    options.algorithm = algo
    return options
end

"""
`options = options_rationals(symmetry, version)`
Initialization of options in rational coordinates. 
"""
function options_rationals(symmetry, version)
    options = options_from_algo(7 + version*8, symmetry, 0, 0) # because in rationals
    options.rational = true
    return options
end

#-----------------------------------------------------------------------
##### Simple tools 
#-----------------------------------------------------------------------

function val_to_char(x)
    c = ' ' # if x ==0
    if x > 0
        c = '1'
    elseif x < 0
        c = '0'
    end
    return c
end

function arr_to_str(stem)
    return join(map(val_to_char, stem))
end

function val_to_char_bis(x)
    c = '0'
    if x > 0
        c ='1'
    end
    return c
end

function arr_to_str_bis(stem)
    return join(map(val_to_char_bis, stem))
end

function comp(value)
    return (value == 1) ? 0 : 1
end

function opposite(vector)
    return map(comp, vector)
end

function symmetric_count(svs)
    n = length(svs)
    p = length(svs[1])
    count = 0
    checked = zeros(Int, n)
    for i in 1:n-1
        if checked[i] == 0
            checked[i] = 1
            truci = svs[i]
            for j in i+1:n
                if truci == -svs[j]
                    count += 2
                    checked[j] = 1
                    break
                end
            end
        end
    end
    return count
end