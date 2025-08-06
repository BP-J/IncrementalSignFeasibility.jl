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

### File used to run some benchmark and comparisons - does not contain any function.

using DelimitedFiles
using Printf

list_instances_rand = ["rand_4_8_2.txt", "rand_7_8_4.txt", "rand_7_9_4.txt", "rand_7_10_5.txt", "rand_7_11_4.txt", "rand_7_12_6.txt", "rand_7_13_5.txt", "rand_7_14_7.txt", "rand_8_15_7.txt", "rand_9_16_8.txt", "rand_10_17_9.txt"]
list_instances_2d = ["degen2d_20_4.txt", "degen2d_20_5.txt", "degen2d_20_6.txt", "degen2d_20_7.txt", "degen2d_20_8.txt"]
list_instances_srand = ["srand_8_20_2.txt", "srand_8_20_4.txt", "srand_8_20_6.txt"]
list_instances_perm = ["perm_5.txt", "perm_6.txt", "perm_7.txt", "perm_8.txt"]
list_instances_ratio = ["ratio_20_3_7.txt", "ratio_20_3_9.txt", "ratio_20_4_7.txt", "ratio_20_4_9.txt", "ratio_20_5_7.txt", "ratio_20_5_9.txt", "ratio_20_6_7.txt", "ratio_20_6_9.txt", "ratio_20_7_7.txt", "ratio_20_7_9.txt"]

list_instances_threshold = ["BEK-threshold_3.txt", "BEK-threshold_4.txt", "BEK-threshold_5.txt"]
list_instances_resonance = ["BEK-resonance_4.txt", "BEK-resonance_5.txt"] #, "BEK-resonance_6.txt"
list_instances_demicube = ["BEK-demicube_4.txt", "BEK-demicube_5.txt", "BEK-demicube_6.txt"]
list_instances_crosspoly = ["BEK-crosspoly_5.txt", "BEK-crosspoly_6.txt", "BEK-crosspoly_7.txt", "BEK-crosspoly_8.txt", "BEK-crosspoly_9.txt", "BEK-crosspoly_10.txt", "BEK-crosspoly_11.txt", "BEK-crosspoly_12.txt"]


options_0 = Options(0, 0, false, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true)
options_1 = Options(1, 0, false, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true)
options_2 = Options(2, 0, true, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true)
options_3 = Options(3, 3, true, false, true, 0, false, true, 100000*eps(), 1000*eps(), false, true)
options_4 = Options(4, 3, true, false, true, 1, false, true, 100000*eps(), 1000*eps(), false, true)
options_4r = Options(4, 3, true, false, true, 1, true, true, 100000*eps(), 1000*eps(), false, true)    # unused, no recursive covering when very few stem vectors 
options_5 = Options(5, 3, true, false, true, 2, false, true, 100000*eps(), 1000*eps(), false, true)
options_5r = Options(5, 3, true, false, true, 2, true, true, 100000*eps(), 1000*eps(), false, true)    # unused, no recursive covering when very few stem vectors 
options_6 = Options(6, 3, true, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, true)    # w often faster
options_6r = Options(6, 3, true, false, true, 3, true, true, 100000*eps(), 1000*eps(), false, true)    # w often faster
options_6w = Options(6, 3, true, false, true, 3, false, true, 100000*eps(), 1000*eps(), true, true)
options_6wr = Options(6, 3, true, false, true, 3, true, true, 100000*eps(), 1000*eps(), true, true)
options_7 = Options(7, 0, false, false, true, 3, false, true, 100000*eps(), 1000*eps(), false, false)  # w often faster
options_7r = Options(7, 0, false, false, true, 3, true, true, 100000*eps(), 1000*eps(), false, false)  # w often faster
options_7w = Options(7, 0, false, false, true, 3, false, true, 100000*eps(), 1000*eps(), true, false)
options_7wr = Options(7, 0, false, false, true, 3, true, true, 100000*eps(), 1000*eps(), true, false)

options_l = Options(7, 0, false, false, true, 4, true, true, 100000*eps(), 1000*eps(), true, false)    # no stree structure

algos = [0,1,2,3,4,5,6,7]
# algos = ["l"]

# list_instances = [list_instances_rand ; list_instances_2d ; list_instances_srand ; list_instances_perm ; list_instances_ratio]
list_instances = ["BEK-resonance_6.txt"]

### end of the comparison for the recursive or not covering shenanigans

for name in list_instances

    if name[1:4] == "BEK-"
        fullname = "linear_data/" * name[5:end]
    else
        fullname = "../../complementarity-internship-2021/papers/1-bdiffmin-affine/matlab/isf/data/data_" * name
    end
    matrice = readdlm(fullname)
    
    println(name)

    if 0 in algos
        info_0 = isf(matrice, options_0)
        time_0 = info_0.cput_total
        time_lp_0 = info_0.cput_lp
        list_times_0 = [time_0]
        list_tlp_0 = [time_lp_0]
        count_0 = 1
        while (count_0 < 2 && time_0 < 1800) || (time_0 < 100 && count_0 < 100)
            count_0 += 1
            info_0 = isf(matrice, options_0)
            time_0 += info_0.cput_total
            time_lp_0 += info_0.cput_lp
            append!(list_times_0, info_0.cput_total)
            append!(list_tlp_0, info_0.cput_lp)
        end
        avg_time_0 = time_0 / count_0
        avg_lp_0 = time_lp_0 / count_0
        #println("Average time for $(count_0) trials, time_0 = $(avg_time_0), time_lp_0 = $(avg_lp_0)")
        # println(list_times_0)
        # println(list_tlp_0)
        string_avg_0 = @sprintf "%.2E" avg_time_0
        string_avg_lp_0 = @sprintf "%.2E" avg_lp_0
        # println("\\texttt{$(name)} & \$" * string_avg_0 * "\$ & \$" * string_avg_lp_0 * "\$ \\\\ \\hline")
    end

    ##################################################
    ##################################################

    if 1 in algos
        info_1 = isf(matrice, options_1)
        time_1 = info_1.cput_total
        time_lp_1 = info_1.cput_lp
        list_times_1 = [time_1]
        list_tlp_1 = [time_lp_1]
        count_1 = 1
        while (count_1 < 2 && time_1 < 1800) || (time_1 < 100 && count_1 < 100)
            count_1 += 1
            info_1 = isf(matrice, options_1)
            time_1 += info_1.cput_total
            time_lp_1 += info_1.cput_lp
            append!(list_times_1, info_1.cput_total)
            append!(list_tlp_1, info_1.cput_lp)
        end
        avg_time_1 = time_1 / count_1
        avg_lp_1 = time_lp_1 / count_1
        #println("Average time for $(count_1) trials, time_1 = $(avg_time_1), time_lp_1 = $(avg_lp_1)")
        # println(list_times_1)
        # println(list_tlp_1)
        string_avg_1 = @sprintf "%.2E" avg_time_1
        string_avg_lp_1 = @sprintf "%.2E" avg_lp_1
        string_ratio_1 = @sprintf "%.2E" avg_time_0 / avg_time_1
        string_ratio_lp_1 = @sprintf "%.2E" avg_lp_1 / avg_time_1
    end

    ##################################################
    ##################################################

    if 2 in algos
        info_2 = isf(matrice, options_2)
        time_2 = info_2.cput_total
        time_lp_2 = info_2.cput_lp
        list_times_2 = [time_2]
        list_tlp_2 = [time_lp_2]
        count_2 = 1
        while (count_2 < 2 && time_2 < 1800) || (time_2 < 100 && count_2 < 100)
            count_2 += 1
            info_2 = isf(matrice, options_2)
            time_2 += info_2.cput_total
            time_lp_2 += info_2.cput_lp
            append!(list_times_2, info_2.cput_total)
            append!(list_tlp_2, info_2.cput_lp)
        end
        avg_time_2 = time_2 / count_2
        avg_lp_2 = time_lp_2 / count_2
        #println("Average time for $(count_2) trials, time_2 = $(avg_time_2), time_lp_2 = $(avg_lp_2)")
        # println(list_times_2)
        # println(list_tlp_2)
        string_avg_2 = @sprintf "%.2E" avg_time_2
        string_avg_lp_2 = @sprintf "%.2E" avg_lp_2
        string_ratio_2 = @sprintf "%.2E" avg_time_0 / avg_time_2
        string_ratio_lp_2 = @sprintf "%.2E" avg_lp_2 / avg_time_2
    end

    ##################################################
    ##################################################

    if 3 in algos
        info_3 = isf(matrice, options_3)
        time_3 = info_3.cput_total
        time_lp_3 = info_3.cput_lp
        list_times_3 = [time_3]
        list_tlp_3 = [time_lp_3]
        count_3 = 1
        while (count_3 < 2 && time_3 < 1800) || (time_3 < 100 && count_3 < 100)
            count_3 += 1
            info_3 = isf(matrice, options_3)
            time_3 += info_3.cput_total
            time_lp_3 += info_3.cput_lp
            append!(list_times_3, info_3.cput_total)
            append!(list_tlp_3, info_3.cput_lp)
        end
        avg_time_3 = time_3 / count_3
        avg_lp_3 = time_lp_3 / count_3
        #println("Average time for $(count_3) trials, time_3 = $(avg_time_3), time_lp_3 = $(avg_lp_3)")
        # println(list_times_3)
        # println(list_tlp_3)
        string_avg_3 = @sprintf "%.2E" avg_time_3
        string_avg_lp_3 = @sprintf "%.2E" avg_lp_3
        string_ratio_3 = @sprintf "%.2E" avg_time_0 / avg_time_3
        string_ratio_lp_3 = @sprintf "%.2E" avg_lp_3 / avg_time_3
        # println("\\texttt{$(name)} & \$" * string_avg_0 * "\$ & \$" * string_avg_lp_0 * "\$ & \$" * string_avg_3 * "\$ & \$" * string_avg_lp_3 * "\$ \\\\ \\hline")
    end
    
    ##################################################
    ##################################################

    if 4 in algos
        info_4 = isf(matrice, options_4)
        time_4 = info_4.cput_total
        time_lp_4 = info_4.cput_lp
        time_sv_4 = info_4.cput_sv
        time_cover_4 = info_4.cput_cover
        list_times_4 = [time_4]
        list_tlp_4 = [time_lp_4]
        list_sv_4 = [time_sv_4]
        list_cover_4 = [time_cover_4]
        count_4 = 1
        while (count_4 < 2 && time_4 < 1800) || (time_4 < 100 && count_4 < 100)
            count_4 += 1
            info_4 = isf(matrice, options_4)
            time_4 += info_4.cput_total
            time_lp_4 += info_4.cput_lp
            time_sv_4 += info_4.cput_sv
            time_cover_4 += info_4.cput_cover
            append!(list_times_4, info_4.cput_total)
            append!(list_tlp_4, info_4.cput_lp)
            append!(list_sv_4, info_4.cput_sv)
            append!(list_cover_4, info_4.cput_cover)
        end
        avg_time_4 = time_4 / count_4
        avg_lp_4 = time_lp_4 / count_4
        avg_sv_4 = time_sv_4 / count_4
        avg_cover_4 = time_cover_4 / count_4
        #println("Average time for $(count_4) trials, time_4 = $(avg_time_4), time_lp_4 = $(avg_lp_4), time_sv_4 = $(avg_sv_4), time_cover_4 = $(avg_cover_4)")
        # println(list_times_4)
        # println(list_tlp_4)
        string_avg_4 = @sprintf "%.2E" avg_time_4
        string_avg_lp_4 = @sprintf "%.2E" avg_lp_4
        string_ratio_4 = @sprintf "%.2E" avg_time_0 / avg_time_4
        string_ratio_lp_4 = @sprintf "%.2E" avg_lp_4 / avg_time_4
        string_sv_4 = @sprintf "%.2E" avg_sv_4
        string_cover_4 = @sprintf "%.2E" avg_cover_4
        string_ratio_sv_4 = @sprintf "%.2E" avg_sv_4 / avg_time_4
        string_ratio_cover_4 = @sprintf "%.2E" avg_cover_4 / avg_time_4
    end
    
    ##################################################
    ##################################################

    if 5 in algos
        info_5 = isf(matrice, options_5)
        time_5 = info_5.cput_total
        time_lp_5 = info_5.cput_lp
        time_sv_5 = info_5.cput_sv
        time_cover_5 = info_5.cput_cover
        list_times_5 = [time_5]
        list_tlp_5 = [time_lp_5]
        list_sv_5 = [time_sv_5]
        list_cover_5 = [time_cover_5]
        count_5 = 1
        while (count_5 < 2 && time_5 < 1800) || (time_5 < 100 && count_5 < 100)
            count_5 += 1
            info_5 = isf(matrice, options_5)
            time_5 += info_5.cput_total
            time_lp_5 += info_5.cput_lp
            time_sv_5 += info_5.cput_sv
            time_cover_5 += info_5.cput_cover
            append!(list_times_5, info_5.cput_total)
            append!(list_tlp_5, info_5.cput_lp)
            append!(list_sv_5, info_5.cput_sv)
            append!(list_cover_5, info_5.cput_cover)
        end
        avg_time_5 = time_5 / count_5
        avg_lp_5 = time_lp_5 / count_5
        avg_sv_5 = time_sv_5 / count_5
        avg_cover_5 = time_cover_5 / count_5
        #println("Average time for $(count_5) trials, time_5 = $(avg_time_5), time_lp_5 = $(avg_lp_5), time_sv_5 = $(avg_sv_5), time_cover_5 = $(avg_cover_5)")
        # println(list_times_5)
        # println(list_tlp_5)
        string_avg_5 = @sprintf "%.2E" avg_time_5
        string_avg_lp_5 = @sprintf "%.2E" avg_lp_5
        string_ratio_5 = @sprintf "%.2E" avg_time_0 / avg_time_5
        string_ratio_lp_5 = @sprintf "%.2E" avg_lp_5 / avg_time_5
        string_sv_5 = @sprintf "%.2E" avg_sv_5
        string_cover_5 = @sprintf "%.2E" avg_cover_5
        string_ratio_sv_5 = @sprintf "%.2E" avg_sv_5 / avg_time_5
        string_ratio_cover_5 = @sprintf "%.2E" avg_cover_5 / avg_time_5
    end
    
    ##################################################
    ##################################################

    if 6 in algos
        # options_for_6 = Options(6, 3, true, 1, 1, true, 0, 3, false, true, 0, 0, true, true)
        # if name in ["rand_10_17_9.txt", "srand_8_20_4.txt", "srand_8_20_6.txt", "perm_8.txt", "ratio_20_4_7.txt", "ratio_20_4_9.txt", "ratio_20_5_7.txt", "ratio_20_5_9.txt", "ratio_20_6_7.txt", "ratio_20_6_9.txt", "ratio_20_7_7.txt", "ratio_20_7_9.txt"]
        #     options_for_6.svsprod = true    # for these ones faster to use the recursive test
        # end
        # if name in ["rand_4_8_2.txt", "degen2d_20_4.txt", "degen2d_20_5.txt", "degen2d_20_6.txt", "degen2d_20_7.txt", "degen2d_20_8.txt", "ratio_20_3_7.txt", "ratio_20_3_9.txt"]
        #     options_for_6.wechelon = false  # for these ones faster NOT to use the echelon form (so by default we take it because there are less instances here)
        # end

        # info_6 = isf(matrice, options_6)
        # time_6 = info_6.cput_total
        # time_lp_6 = info_6.cput_lp
        # time_sv_6 = info_6.cput_sv
        # time_cover_6 = info_6.cput_cover
        # list_times_6 = [time_6]
        # list_tlp_6 = [time_lp_6]
        # list_sv_6 = [time_sv_6]
        # list_cover_6 = [time_cover_6]
        # count_6 = 1
        # while (count_6 < 2 && time_6 < 1800) || (time_6 < 100 && count_6 < 100)
        #     count_6 += 1
        #     info_6 = isf(matrice, options_6)
        #     time_6 += info_6.cput_total
        #     time_lp_6 += info_6.cput_lp
        #     time_sv_6 += info_6.cput_sv
        #     time_cover_6 += info_6.cput_cover
        #     append!(list_times_6, info_6.cput_total)
        #     append!(list_tlp_6, info_6.cput_lp)
        #     append!(list_sv_6, info_6.cput_sv)
        #     append!(list_cover_6, info_6.cput_cover)
        # end
        # avg_time_6 = time_6 / count_6
        # avg_lp_6 = time_lp_6 / count_6
        # avg_sv_6 = time_sv_6 / count_6
        # avg_cover_6 = time_cover_6 / count_6
        # #println("Average time for $(count_6) trials, time_6 = $(avg_time_6), time_lp_6 = $(avg_lp_6), time_sv_6 = $(avg_sv_6), time_cover_6 = $(avg_cover_6)")
        # # println(list_times_6)
        # # println(list_tlp_6)
        # string_avg_6 = @sprintf "%.2E" avg_time_6
        # string_avg_lp_6 = @sprintf "%.2E" avg_lp_6
        # string_ratio_6 = @sprintf "%.2E" avg_time_0 / avg_time_6
        # string_ratio_lp_6 = @sprintf "%.2E" avg_lp_6 / avg_time_6
        # string_sv_6 = @sprintf "%.2E" avg_sv_6
        # string_cover_6 = @sprintf "%.2E" avg_cover_6
        # string_ratio_sv_6 = @sprintf "%.2E" avg_sv_6 / avg_time_6
        # string_ratio_cover_6 = @sprintf "%.2E" avg_cover_6 / avg_time_6
        string_avg_6 = @sprintf "*" 
        string_avg_lp_6 = @sprintf "*" 
        string_ratio_6 = @sprintf "*" 
        string_ratio_lp_6 = @sprintf "*" 
        string_sv_6 = @sprintf "*" 
        string_cover_6 = @sprintf "*" 
        string_ratio_sv_6 = @sprintf "*" 
        string_ratio_cover_6 = @sprintf "*" 
    end

    ##################################################
    ##################################################

    if 7 in algos 
        # options_for_7 = Options(7, 0, false, 1, 1, true, 0, 3, true, true, 0, 0, true, false)
        # if name in ["rand_9_16_8.txt", "rand_10_17_9.txt", "srand_8_20_4.txt", "srand_8_20_6.txt", "perm_7.txt", "perm_8.txt", "ratio_20_4_7.txt", "ratio_20_4_9.txt", "ratio_20_5_7.txt", "ratio_20_5_9.txt", "ratio_20_6_7.txt", "ratio_20_6_9.txt", "ratio_20_7_7.txt", "ratio_20_7_9.txt"]
        #     options_for_7.svsprod = true
        # end
        # if name in ["rand_4_8_2.txt", "degen2d_20_4.txt", "degen2d_20_5.txt", "degen2d_20_6.txt", "degen2d_20_7.txt", "degen2d_20_8.txt", "ratio_20_3_7.txt", "ratio_20_3_9.txt"]
        #     options_for_7.wechelon = false  # for these ones faster NOT to use the echelon form (so by default we take it because there are less instances here)
        # end

        # info_7 = isf(matrice, options_7)
        # time_7 = info_7.cput_total
        # time_lp_7 = info_7.cput_lp
        # time_sv_7 = info_7.cput_sv
        # time_cover_7 = info_7.cput_cover
        # list_times_7 = [time_7]
        # list_tlp_7 = [time_lp_7]
        # list_sv_7 = [time_sv_7]
        # list_cover_7 = [time_cover_7]
        # count_7 = 1
        # while (count_7 < 2 && time_7 < 1800) || (time_7 < 100 && count_7 < 100)
        #     count_7 += 1
        #     info_7 = isf(matrice, options_7)
        #     time_7 += info_7.cput_total
        #     time_lp_7 += info_7.cput_lp
        #     time_sv_7 += info_7.cput_sv
        #     time_cover_7 += info_7.cput_cover
        #     append!(list_times_7, info_7.cput_total)
        #     append!(list_tlp_7, info_7.cput_lp)
        #     append!(list_sv_7, info_7.cput_sv)
        #     append!(list_cover_7, info_7.cput_cover)
        # end
        # avg_time_7 = time_7 / count_7
        # avg_lp_7 = time_lp_7 / count_7
        # avg_sv_7 = time_sv_7 / count_7
        # avg_cover_7 = time_cover_7 / count_7
        # #println("Average time for $(count_7) trials, time_7 = $(avg_time_7), time_lp_7 = $(avg_lp_7), time_sv_7 = $(avg_sv_7), time_cover_7 = $(avg_cover_7)")
        # # println(list_times_7)
        # # println(list_tlp_7)
        # string_avg_7 = @sprintf "%.2E" avg_time_7
        # string_avg_lp_7 = @sprintf "%.2E" avg_lp_7                  # 0 but whatever
        # string_ratio_7 = @sprintf "%.2E" avg_time_0 / avg_time_7
        # string_ratio_lp_7 = @sprintf "%.2E" avg_lp_7 / avg_time_7   # 0 but whatever
        # string_sv_7 = @sprintf "%.2E" avg_sv_7
        # string_cover_7 = @sprintf "%.2E" avg_cover_7
        # string_ratio_sv_7 = @sprintf "%.2E" avg_sv_7 / avg_time_7
        # string_ratio_cover_7 = @sprintf "%.2E" avg_cover_7 / avg_time_7
        # println("\\texttt{$(name)} & \$" * string_avg_7 * "\$ & \$" * string_cover_7 * "\$ \\\\ \\hline")
        string_avg_7 = @sprintf "*" 
        string_avg_lp_7 = @sprintf "*" 
        string_ratio_7 = @sprintf "*" 
        string_ratio_lp_7 = @sprintf "*" 
        string_sv_7 = @sprintf "*" 
        string_cover_7 = @sprintf "*" 
        string_ratio_sv_7 = @sprintf "*"
        string_ratio_cover_7 = @sprintf "*" 
    end

    ##################################################
    ##################################################

    if "l" in algos
        info_l = isf(matrice, options_l)
        time_l = info_l.cput_total
        time_lp_l = info_l.cput_lp
        time_sv_l = info_l.cput_sv
        time_cover_l = info_l.cput_cover
        list_times_l = [time_l]
        list_tlp_l = [time_lp_l]
        list_sv_l = [time_sv_l]
        list_cover_l = [time_cover_l]
        count_l = 1
        while (count_l < 2 && time_l < 1800) || (time_l < 100 && count_l < 100)
            count_l += 1
            info_l = isf(matrice, options_l)
            time_l += info_l.cput_total
            time_lp_l += info_l.cput_lp
            time_sv_l += info_l.cput_sv
            time_cover_l += info_l.cput_cover
            append!(list_times_l, info_l.cput_total)
            append!(list_tlp_l, info_l.cput_lp)
            append!(list_sv_l, info_l.cput_sv)
            append!(list_cover_l, info_l.cput_cover)
        end
        avg_time_l = time_l / count_l
        avg_lp_l = time_lp_l / count_l
        avg_sv_l = time_sv_l / count_l
        avg_cover_l = time_cover_l / count_l
        string_avg_l = @sprintf "%.2E" avg_time_l
        string_avg_lp_l = @sprintf "%.2E" avg_lp_l
        string_sv_l = @sprintf "%.2E" avg_sv_l
        string_cover_l = @sprintf "%.2E" avg_cover_l

        println("$(name)")
        println("time = " * string_avg_l * ", time sv = " * string_sv_l * ", time cover = " * string_cover_l * ", time lp = " * string_avg_lp_l)
    end


    ##################################################


    if algos == [0,1,2,3,4,5,6,7]
        println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_0 * "\$} & \$" * string_avg_lp_0 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_1 * "\$} & \$" * string_avg_lp_1 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_2 * "\$} & \$" * string_avg_lp_2 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_3 * "\$} & \$" * string_avg_lp_3 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_4 * "\$} & \$" * string_avg_lp_4 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_5 * "\$} & \$" * string_avg_lp_5 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_6 * "\$} & \$" * string_avg_lp_6 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_7 * "\$} & \$" * string_avg_lp_7 * "\$ \\\\ \\hline")
        println("ratio | \\% & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{\$" * string_ratio_1 * "\$} & \$" * string_ratio_lp_1 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_2 * "\$} & \$" * string_ratio_lp_2 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_3 * "\$} & \$" * string_ratio_lp_3 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_4 * "\$} & \$" * string_ratio_lp_4 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_5 * "\$} & \$" * string_ratio_lp_5 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_6 * "\$} & \$" * string_ratio_lp_6 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_7 * "\$} & \$" * string_ratio_lp_7 * "\$ \\\\ \\hline")
        println("sv | cover & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{\$" * string_sv_4 * "\$} & \$" * string_cover_4 * "\$ & \\multicolumn{1}{c|}{\$" * string_sv_5 * "\$} & \$" * string_cover_5 * "\$ & \\multicolumn{1}{c|}{\$" * string_sv_6 * "\$} & \$" * string_cover_6 * "\$ & \\multicolumn{1}{c|}{\$" * string_sv_7 * "\$} & \$" * string_cover_7 * "\$ \\\\ \\hline")
        println("\\% | \\% & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{\$" * string_ratio_sv_4 * "\$} & \$" * string_ratio_cover_4 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_sv_5 * "\$} & \$" * string_ratio_cover_5 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_sv_6 * "\$} & \$" * string_ratio_cover_6 * "\$ & \\multicolumn{1}{c|}{\$" * string_ratio_sv_7 * "\$} & \$" * string_ratio_cover_7 * "\$ \\\\ \\hline\\hline")
        println()
        println("autre ligne")
        println("Problem & \\multicolumn{2}{c|}{RC} & \\multicolumn{2}{c|}{A} & \\multicolumn{2}{c|}{AB} & \\multicolumn{2}{c|}{ABC} & \\multicolumn{2}{c|}{D1} & \\multicolumn{2}{c|}{D2} & \\multicolumn{2}{c|}{D3} & \\multicolumn{2}{c|}{D4} \\\\ \\hline")
        println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_0 * "\$} & " * " " * " & \\multicolumn{1}{c|}{\$" * string_avg_1 * "\$} & {\\color{blue}{" * string_ratio_1* "}}& \\multicolumn{1}{c|}{\$" * string_avg_2 * "\$} & {\\color{blue}{" * string_ratio_2* "}} & \\multicolumn{1}{c|}{\$" * string_avg_3 * "\$} & {\\color{blue}{" * string_ratio_3* "}}& \\multicolumn{1}{c|}{\$" * string_avg_4 * "\$} & {\\color{blue}{" * string_ratio_4* "}} & \\multicolumn{1}{c|}{\$" * string_avg_5 * "\$} & {\\color{blue}{" * string_ratio_5* "}}& \\multicolumn{1}{c|}{\$" * string_avg_6 * "\$} & {\\color{blue}{" * string_ratio_6* "}} & \\multicolumn{1}{c|}{\$" * string_avg_7 * "\$} & {\\color{blue}{" * string_ratio_7* "}}\\\\ \\hline")
    end
end

