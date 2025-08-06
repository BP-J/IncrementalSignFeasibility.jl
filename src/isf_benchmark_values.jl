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
`... = isf_benchmark_values(info, options)`
Function returning some information about some instance values, see `isf_benchmark_launches_aff.jl`. 
"""
function isf_benchmark_values(info::Info, options::Options)

    time = info.cput_total

    FLP = info.nb_feaslop
    ILP = info.nb_infeaslop
    time_LP = info.cput_lp    
    if options.algorithm in [7,15]
        avg_LP = 0
    else
        avg_LP = time_LP / (FLP + ILP)
    end
    prop_LP = time_LP / time

    syms = info.nb_stems_sym
    asyms = info.nb_stems_asym
    dupli_syms = info.nb_duplicated_stems_sym
    dupli_asyms = info.nb_duplicated_stems_asym
    time_sv = info.cput_sv
    prop_sv = time_sv / time

    checks = info.nb_sv_checks
    time_cover = info.cput_cover
    if mod(options.algorithm, 8) <= 3
        detect_ratio = 0
        avg_cover = 0
    else
        detect_ratio = info.nb_sv_detect / checks
        avg_cover = time_cover / checks
    end
    prop_cover = time_cover / time

    return time, FLP, ILP, time_LP, avg_LP, prop_LP, syms, asyms, dupli_syms, dupli_asyms, time_sv, prop_sv, checks, detect_ratio, time_cover, avg_cover, prop_cover
end