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

list_instances_rand = ["R_8_2", "R_8_4", "R_9_4", "R_10_5", "R_11_4", "R_12_6", "R_13_5", "R_14_7", "R_15_7", "R_16_8", "R_17_9"]
list_instances_2d = ["degen2d_4", "degen2d_5", "degen2d_6", "degen2d_7", "degen2d_8"]
list_instances_srand = ["srand_8_20_2", "srand_8_20_4", "srand_8_20_6"]
list_instances_perm = ["perm_5", "perm_6", "perm_7", "perm_8"]
list_instances_ratio = ["ratio_3_7", "ratio_3_9", "ratio_4_7", "ratio_4_9", "ratio_5_7", "ratio_5_9", "ratio_6_7", "ratio_6_9", "ratio_7_7", "ratio_7_9"]

list_instances = [list_instances_rand ; list_instances_2d ; list_instances_srand ; list_instances_perm ; list_instances_ratio]

function benchmark_affine(list_instances)

    options_0 = Options(0, 0, false, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_1 = Options(1, 0, false, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_2 = Options(2, 0, true, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_3 = Options(3, 3, true, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_4 = Options(4, 3, true, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, true)
    options_4r = Options(4, 3, true, false, true, 1, true, false, 100000*eps(), 1000*eps(), false, true)    # unused, no recursive covering when very few stem vectors 
    options_5 = Options(5, 3, true, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, true)
    options_5r = Options(5, 3, true, false, true, 2, true, false, 100000*eps(), 1000*eps(), false, true)    # unused, no recursive covering when very few stem vectors 
    options_6 = Options(6, 3, true, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, true)    # w often faster
    options_6r = Options(6, 3, true, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, true)    # w often faster
    options_6w = Options(6, 3, true, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, true)
    options_6wr = Options(6, 3, true, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, true)
    options_7 = Options(7, 0, false, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false)  # w often faster
    options_7r = Options(7, 0, false, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false)  # w often faster
    options_7w = Options(7, 0, false, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false)
    options_7wr = Options(7, 0, false, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false)

    options_8 = Options(8, 0, false, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_9 = Options(9, 0, false, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_10 = Options(10, 0, true, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_11 = Options(11, 3, true, false, true, 0, false, false, 100000*eps(), 1000*eps(), false, true)
    options_12 = Options(12, 3, true, false, true, 1, false, false, 100000*eps(), 1000*eps(), false, true)
    options_12r = Options(12, 3, true, false, true, 1, true, false, 100000*eps(), 1000*eps(), false, true)  # unused, no recursive covering when very few stem vectors 
    options_13 = Options(13, 3, true, false, true, 2, false, false, 100000*eps(), 1000*eps(), false, true)
    options_13r = Options(13, 3, true, false, true, 2, true, false, 100000*eps(), 1000*eps(), false, true)  # unused, no recursive covering when very few stem vectors 
    options_14 = Options(14, 3, true, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, true)  # w often faster
    options_14r = Options(14, 3, true, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, true)  # w often faster
    options_14w = Options(14, 3, true, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, true)
    options_14wr = Options(14, 3, true, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, true)
    options_15 = Options(15, 0, false, false, true, 3, false, false, 100000*eps(), 1000*eps(), false, false)# w often faster
    options_15r = Options(15, 0, false, false, true, 3, true, false, 100000*eps(), 1000*eps(), false, false)# w often faster
    options_15w = Options(15, 0, false, false, true, 3, false, false, 100000*eps(), 1000*eps(), true, false)
    options_15wr = Options(15, 0, false, false, true, 3, true, false, 100000*eps(), 1000*eps(), true, false)

    options_l = Options(7, 0, false, false, true, 4, true, false, 100000*eps(), 1000*eps(), true, false)    # no tree structure

    algos = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

    DataFrame_instance = DataFrame(total_time = Float64[], FeasLP = Int64[], InfeasLP = Int64[], ratio_LP = Float64[], LP_time = Float64[], avg_LP = Float64[], prop_time_LP = Float64[],
                                    syms = Int64[], asyms = Int64[], dupli_syms = Int64[], dupli_asyms = Int64[], sv_time = Float64[], prop_time_sv = Float64[],
                                    checks = Int64[], detect_ratio = Float64[], cover_time = Float64[], avg_cover = Float64[], prop_time_cover = Float64[])

    DataFrames_vector = [DataFrame_instance for i in 1:length(list_instances)]
    # DataFrames_vector_HnH = [DataFrame_instance for i in 1:length(list_instances)]

    DataFrames_times = DataFrame(Name = String[], RC = Float64[], RC_HnH = Float64[], RC_ratio_HnH = Float64[],
                                    A = Float64[], A_ratio = Float64[], A_HnH = Float64[], A_ratio_HnH = Float64[], 
                                    B = Float64[], B_ratio = Float64[], B_HnH = Float64[], B_ratio_HnH = Float64[], 
                                    C = Float64[], C_ratio = Float64[], C_HnH = Float64[], C_ratio_HnH = Float64[], 
                                    D1 = Float64[], D1_ratio = Float64[], D1_HnH = Float64[], D1_ratio_HnH = Float64[], 
                                    D2 = Float64[], D2_ratio = Float64[], D2_HnH = Float64[], D2_ratio_HnH = Float64[], 
                                    D3 = Float64[], D3_ratio = Float64[], D3_HnH = Float64[], D3_ratio_HnH = Float64[], 
                                    D4 = Float64[], D4_ratio = Float64[], D4_HnH = Float64[], D4_ratio_HnH = Float64[])

    for DF_index in eachindex(list_instances)

        name = list_instances[DF_index]
        fullname = "affine_data/data_aff_" * name * ".txt"
        matrice = readdlm(fullname)
        println(name,"\n")

        if 0 in algos
            info_0 = isf(matrice, options_0)
            count_0 = 1

            time_0, FLP_0, ILP_0, time_LP_0, avg_LP_0, prop_LP_0, syms_0, asyms_0, dupli_syms_0, dupli_asyms_0, time_sv_0, prop_sv_0, checks_0, detect_ratio_0, time_cover_0, avg_cover_0, prop_cover_0 = isf_benchmark_values(info_0, options_0)

            while count_0 < 3 || (time_0 < 100 && count_0 < 100)
                count_0 += 1
                info_0 = isf(matrice, options_0)
                time_0       += info_0.cput_total
                time_LP_0    += info_0.cput_lp
                prop_LP_0    += (info_0.cput_lp / info_0.cput_total)
                time_sv_0    += info_0.cput_sv
                prop_sv_0    += (info_0.cput_sv / info_0.cput_total)
                time_cover_0 += info_0.cput_cover
                prop_cover_0 += (info_0.cput_cover / info_0.cput_total)
            end
            time_0       /= count_0
            time_LP_0    /= count_0
            prop_LP_0    /= count_0
            time_sv_0    /= count_0
            prop_sv_0    /= count_0
            time_cover_0 /= count_0
            prop_cover_0 /= count_0

            push!(DataFrames_vector[DF_index], [time_0, FLP_0, ILP_0, 1.0, time_LP_0, time_LP_0 / (FLP_0 + ILP_0), prop_LP_0, syms_0, asyms_0, dupli_syms_0, dupli_asyms_0, time_sv_0, prop_sv_0, checks_0, detect_ratio_0, time_cover_0, 0, prop_cover_0])
            # println("Average time for $(count_0) trials, time_0 = $(avg_time_0), time_LP_0 = $(avg_LP_0)")
            # string_avg_0 = @sprintf "%.2E" avg_time_0
            # string_avg_LP_0 = @sprintf "%.2E" avg_LP_0
            # string_ppt_LP_0 = @sprintf "%.2f" avg_LP_0 / avg_time_0
        end

        ##################################################
        ##################################################

        if 1 in algos
            info_1 = isf(matrice, options_1)
            count_1 = 1

            time_1, FLP_1, ILP_1, time_LP_1, avg_LP_1, prop_LP_1, syms_1, asyms_1, dupli_syms_1, dupli_asyms_1, time_sv_1, prop_sv_1, checks_1, detect_ratio_1, time_cover_1, avg_cover_1, prop_cover_1 = isf_benchmark_values(info_1, options_1)

            while count_1 < 3 || (time_1 < 100 && count_1 < 100)
                count_1 += 1
                info_1 = isf(matrice, options_1)
                time_1       += info_1.cput_total
                time_LP_1    += info_1.cput_lp
                prop_LP_1    += (info_1.cput_lp / info_1.cput_total)
                time_sv_1    += info_1.cput_sv
                prop_sv_1    += (info_1.cput_sv / info_1.cput_total)
                time_cover_1 += info_1.cput_cover
                prop_cover_1 += (info_1.cput_cover / info_1.cput_total)
            end
            time_1       /= count_1
            time_LP_1    /= count_1
            prop_LP_1    /= count_1
            time_sv_1    /= count_1
            prop_sv_1    /= count_1
            time_cover_1 /= count_1
            prop_cover_1 /= count_1

            push!(DataFrames_vector[DF_index], [time_1, FLP_1, ILP_1, (FLP_0+ILP_0)/(FLP_1+ILP_1) , time_LP_1, time_LP_1 / (FLP_1 + ILP_1), prop_LP_1, syms_1, asyms_1, dupli_syms_1, dupli_asyms_1, time_sv_1, prop_sv_1, checks_1, detect_ratio_1, time_cover_1, 0, prop_cover_1])
            # println("Average time for $(count_1) trials, time_1 = $(avg_time_1), time_LP_1 = $(avg_LP_1)")
            # string_avg_1 = @sprintf "%.2E" avg_time_1
            # string_avg_LP_1 = @sprintf "%.2E" avg_LP_1
            # string_ppt_LP_1 = @sprintf "%.2f" avg_LP_1 / avg_time_1
            # string_ratio_1 = @sprintf "%.2f" avg_time_0 / avg_time_1
            # string_ratio_LP_1 = @sprintf "%.2f" avg_LP_0 / avg_LP_1
            # string_ratio_ppt_LP_1 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_1 / avg_time_1)
        end

        ##################################################
        ##################################################

        if 2 in algos
            info_2 = isf(matrice, options_2)
            count_2 = 1
            
            time_2, FLP_2, ILP_2, time_LP_2, avg_LP_2, prop_LP_2, syms_2, asyms_2, dupli_syms_2, dupli_asyms_2, time_sv_2, prop_sv_2, checks_2, detect_ratio_2, time_cover_2, avg_cover_2, prop_cover_2 = isf_benchmark_values(info_2, options_2)
            
            while count_2 < 3 || (time_2 < 100 && count_2 < 100)
                count_2 += 1
                info_2 = isf(matrice, options_2)
                time_2       += info_2.cput_total
                time_LP_2    += info_2.cput_lp
                prop_LP_2    += (info_2.cput_lp / info_2.cput_total)
                time_sv_2    += info_2.cput_sv
                prop_sv_2    += (info_2.cput_sv / info_2.cput_total)
                time_cover_2 += info_2.cput_cover
                prop_cover_2 += (info_2.cput_cover / info_2.cput_total)
            end
            time_2       /= count_2
            time_LP_2    /= count_2
            prop_LP_2    /= count_2
            time_sv_2    /= count_2
            prop_sv_2    /= count_2
            time_cover_2 /= count_2
            prop_cover_2 /= count_2
            
            push!(DataFrames_vector[DF_index], [time_2, FLP_2, ILP_2, (FLP_0+ILP_0)/(FLP_2+ILP_2) , time_LP_2, time_LP_2 / (FLP_2 + ILP_2), prop_LP_2, syms_2, asyms_2, dupli_syms_2, dupli_asyms_2, time_sv_2, prop_sv_2, checks_2, detect_ratio_2, time_cover_2, 0, prop_cover_2])
            # println("Average time for $(count_2) trials, time_2 = $(avg_time_2), time_LP_2 = $(avg_LP_2)")
            # string_avg_2 = @sprintf "%.2E" avg_time_2
            # string_avg_LP_2 = @sprintf "%.2E" avg_LP_2
            # string_ppt_LP_2 = @sprintf "%.2f" avg_LP_2 / avg_time_2
            # string_ratio_2 = @sprintf "%.2f" avg_time_0 / avg_time_2
            # string_ratio_LP_2 = @sprintf "%.2f" avg_LP_0 / avg_LP_2
            # string_ratio_ppt_LP_2 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_2 / avg_time_2)
        end

        ##################################################
        ##################################################

        if 3 in algos
            info_3 = isf(matrice, options_3)
            count_3 = 1
            
            time_3, FLP_3, ILP_3, time_LP_3, avg_LP_3, prop_LP_3, syms_3, asyms_3, dupli_syms_3, dupli_asyms_3, time_sv_3, prop_sv_3, checks_3, detect_ratio_3, time_cover_3, avg_cover_3, prop_cover_3 = isf_benchmark_values(info_3, options_3)
            
            while count_3 < 3 || (time_3 < 100 && count_3 < 100)
                count_3 += 1
                info_3 = isf(matrice, options_3)
                time_3       += info_3.cput_total
                time_LP_3    += info_3.cput_lp
                prop_LP_3    += (info_3.cput_lp / info_3.cput_total)
                time_sv_3    += info_3.cput_sv
                prop_sv_3    += (info_3.cput_sv / info_3.cput_total)
                time_cover_3 += info_3.cput_cover
                prop_cover_3 += (info_3.cput_cover / info_3.cput_total)
            end
            time_3       /= count_3
            time_LP_3    /= count_3
            prop_LP_3    /= count_3
            time_sv_3    /= count_3
            prop_sv_3    /= count_3
            time_cover_3 /= count_3
            prop_cover_3 /= count_3
            
            push!(DataFrames_vector[DF_index], [time_3, FLP_3, ILP_3, (FLP_0+ILP_0)/(FLP_3+ILP_3) , time_LP_3, time_LP_3 / (FLP_3 + ILP_3), prop_LP_3, syms_3, asyms_3, dupli_syms_3, dupli_asyms_3, time_sv_3, prop_sv_3, checks_3, detect_ratio_3, time_cover_3, 0, prop_cover_3])
            # println("Average time for $(count_3) trials, time_3 = $(avg_time_3), time_LP_3 = $(avg_LP_3)")
            # string_avg_3 = @sprintf "%.2E" avg_time_3
            # string_avg_LP_3 = @sprintf "%.2E" avg_LP_3
            # string_ppt_LP_3 = @sprintf "%.2f" avg_LP_3 / avg_time_3
            # string_ratio_3 = @sprintf "%.2f" avg_time_0 / avg_time_3
            # string_ratio_LP_3 = @sprintf "%.2f" avg_LP_0 / avg_LP_3
            # string_ratio_ppt_LP_3 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_3 / avg_time_3)
        end

        ##################################################
        ##################################################

        if 4 in algos
            info_4 = isf(matrice, options_4)
            count_4 = 1
            
            time_4, FLP_4, ILP_4, time_LP_4, avg_LP_4, prop_LP_4, syms_4, asyms_4, dupli_syms_4, dupli_asyms_4, time_sv_4, prop_sv_4, checks_4, detect_ratio_4, time_cover_4, avg_cover_4, prop_cover_4 = isf_benchmark_values(info_4, options_4)
            
            while count_4 < 3 || (time_4 < 100 && count_4 < 100)
                count_4 += 1
                info_4 = isf(matrice, options_4)
                time_4       += info_4.cput_total
                time_LP_4    += info_4.cput_lp
                prop_LP_4    += (info_4.cput_lp / info_4.cput_total)
                time_sv_4    += info_4.cput_sv
                prop_sv_4    += (info_4.cput_sv / info_4.cput_total)
                time_cover_4 += info_4.cput_cover
                prop_cover_4 += (info_4.cput_cover / info_4.cput_total)
            end
            time_4       /= count_4
            time_LP_4    /= count_4
            prop_LP_4    /= count_4
            time_sv_4    /= count_4
            prop_sv_4    /= count_4
            time_cover_4 /= count_4
            prop_cover_4 /= count_4
            
            push!(DataFrames_vector[DF_index], [time_4, FLP_4, ILP_4, (FLP_0+ILP_0)/(FLP_4+ILP_4) , time_LP_4, time_LP_4 / (FLP_4 + ILP_4), prop_LP_4, syms_4, asyms_4, dupli_syms_4, dupli_asyms_4, time_sv_4, prop_sv_4, checks_4, detect_ratio_4, time_cover_4, time_cover_4 / checks_4, prop_cover_4])
            # println("Average time for $(count_4) trials, time_4 = $(avg_time_4), time_LP_4 = $(avg_LP_4)")
            # string_avg_4 = @sprintf "%.2E" avg_time_4
            # string_avg_LP_4 = @sprintf "%.2E" avg_LP_4
            # string_ppt_LP_4 = @sprintf "%.2E" avg_LP_4 / avg_time_4
            # string_sv_4 = @sprintf "%.2E" avg_sv_4
            # string_cv_4 = @sprintf "%.2E" avg_cv_4
            # string_ppt_sv_4 = @sprintf "%.2f" avg_sv_4 / avg_time_4
            # string_ppt_cv_4 = @sprintf "%.2f" avg_cv_4 / avg_time_4
            # string_ratio_4 = @sprintf "%.2f" avg_time_0 / avg_time_4
            # string_ratio_LP_4 = @sprintf "%.2f" avg_LP_0 / avg_LP_4
            # string_ratio_ppt_LP_4 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_4 / avg_time_4)
        end
        
        ##################################################
        ##################################################

        if 5 in algos
            info_5 = isf(matrice, options_5)
            count_5 = 1
            
            time_5, FLP_5, ILP_5, time_LP_5, avg_LP_5, prop_LP_5, syms_5, asyms_5, dupli_syms_5, dupli_asyms_5, time_sv_5, prop_sv_5, checks_5, detect_ratio_5, time_cover_5, avg_cover_5, prop_cover_5 = isf_benchmark_values(info_5, options_5)
            
            while count_5 < 3 || (time_5 < 100 && count_5 < 100)
                count_5 += 1
                info_5 = isf(matrice, options_5)
                time_5       += info_5.cput_total
                time_LP_5    += info_5.cput_lp
                prop_LP_5    += (info_5.cput_lp / info_5.cput_total)
                time_sv_5    += info_5.cput_sv
                prop_sv_5    += (info_5.cput_sv / info_5.cput_total)
                time_cover_5 += info_5.cput_cover
                prop_cover_5 += (info_5.cput_cover / info_5.cput_total)
            end
            time_5       /= count_5
            time_LP_5    /= count_5
            prop_LP_5    /= count_5
            time_sv_5    /= count_5
            prop_sv_5    /= count_5
            time_cover_5 /= count_5
            prop_cover_5 /= count_5
            
            push!(DataFrames_vector[DF_index], [time_5, FLP_5, ILP_5, (FLP_0+ILP_0)/(FLP_5+ILP_5) , time_LP_5, time_LP_5 / (FLP_5 + ILP_5), prop_LP_5, syms_5, asyms_5, dupli_syms_5, dupli_asyms_5, time_sv_5, prop_sv_5, checks_5, detect_ratio_5, time_cover_5, time_cover_5 / checks_5, prop_cover_5])
            # println("Average time for $(count_5) trials, time_5 = $(avg_time_5), time_LP_5 = $(avg_LP_5)")
            # string_avg_5 = @sprintf "%.2E" avg_time_5
            # string_avg_LP_5 = @sprintf "%.2E" avg_LP_5
            # string_ppt_LP_5 = @sprintf "%.2E" avg_LP_5 / avg_time_5
            # string_sv_5 = @sprintf "%.2E" avg_sv_5
            # string_cv_5 = @sprintf "%.2E" avg_cv_5
            # string_ppt_sv_5 = @sprintf "%.2f" avg_sv_5 / avg_time_5
            # string_ppt_cv_5 = @sprintf "%.2f" avg_cv_5 / avg_time_5
            # string_ratio_5 = @sprintf "%.2f" avg_time_0 / avg_time_5
            # string_ratio_LP_5 = @sprintf "%.2f" avg_LP_0 / avg_LP_5
            # string_ratio_ppt_LP_5 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_5 / avg_time_5)
        end

        ##################################################
        ##################################################

        if 6 in algos
            info_6 = isf(matrice, options_6)
            count_6 = 1
            
            time_6, FLP_6, ILP_6, time_LP_6, avg_LP_6, prop_LP_6, syms_6, asyms_6, dupli_syms_6, dupli_asyms_6, time_sv_6, prop_sv_6, checks_6, detect_ratio_6, time_cover_6, avg_cover_6, prop_cover_6 = isf_benchmark_values(info_6, options_6)
            
            while count_6 < 3 || (time_6 < 100 && count_6 < 100)
                count_6 += 1
                info_6 = isf(matrice, options_6)
                time_6       += info_6.cput_total
                time_LP_6    += info_6.cput_lp
                prop_LP_6    += (info_6.cput_lp / info_6.cput_total)
                time_sv_6    += info_6.cput_sv
                prop_sv_6    += (info_6.cput_sv / info_6.cput_total)
                time_cover_6 += info_6.cput_cover
                prop_cover_6 += (info_6.cput_cover / info_6.cput_total)
            end
            time_6       /= count_6
            time_LP_6    /= count_6
            prop_LP_6    /= count_6
            time_sv_6    /= count_6
            prop_sv_6    /= count_6
            time_cover_6 /= count_6
            prop_cover_6 /= count_6
            
            push!(DataFrames_vector[DF_index], [time_6, FLP_6, ILP_6, (FLP_0+ILP_0)/(FLP_6+ILP_6) , time_LP_6, time_LP_6 / (FLP_6 + ILP_6), prop_LP_6, syms_6, asyms_6, dupli_syms_6, dupli_asyms_6, time_sv_6, prop_sv_6, checks_6, detect_ratio_6, time_cover_6, time_cover_6 / checks_6, prop_cover_6])
            # println("Average time for $(count_6) trials, time_6 = $(avg_time_6), time_LP_6 = $(avg_LP_6)")
            # string_avg_6 = @sprintf "%.2E" avg_time_6
            # string_avg_LP_6 = @sprintf "%.2E" avg_LP_6
            # string_ppt_LP_6 = @sprintf "%.2E" avg_LP_6 / avg_time_6
            # string_sv_6 = @sprintf "%.2E" avg_sv_6
            # string_cv_6 = @sprintf "%.2E" avg_cv_6
            # string_ppt_sv_6 = @sprintf "%.2f" avg_sv_6 / avg_time_6
            # string_ppt_cv_6 = @sprintf "%.2f" avg_cv_6 / avg_time_6
            # string_ratio_6 = @sprintf "%.2f" avg_time_0 / avg_time_6
            # string_ratio_LP_6 = @sprintf "%.2f" avg_LP_0 / avg_LP_6
            # string_ratio_ppt_LP_6 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_6 / avg_time_6)
        end

        ##################################################
        ##################################################

        if 7 in algos
            info_7 = isf(matrice, options_7)
            count_7 = 1
            
            time_7, FLP_7, ILP_7, time_LP_7, avg_LP_7, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, avg_cover_7, prop_cover_7 = isf_benchmark_values(info_7, options_7)
            
            while count_7 < 3 || (time_7 < 100 && count_7 < 100)
                count_7 += 1
                info_7 = isf(matrice, options_7)
                time_7       += info_7.cput_total
                time_LP_7    += info_7.cput_lp
                prop_LP_7    += (info_7.cput_lp / info_7.cput_total)
                time_sv_7    += info_7.cput_sv
                prop_sv_7    += (info_7.cput_sv / info_7.cput_total)
                time_cover_7 += info_7.cput_cover
                prop_cover_7 += (info_7.cput_cover / info_7.cput_total)
            end
            time_7       /= count_7
            time_LP_7    /= count_7
            prop_LP_7    /= count_7
            time_sv_7    /= count_7
            prop_sv_7    /= count_7
            time_cover_7 /= count_7
            prop_cover_7 /= count_7
            
            push!(DataFrames_vector[DF_index], [time_7, FLP_7, ILP_7, 0, time_LP_7, 0, prop_LP_7, syms_7, asyms_7, dupli_syms_7, dupli_asyms_7, time_sv_7, prop_sv_7, checks_7, detect_ratio_7, time_cover_7, time_cover_7 / checks_7, prop_cover_7])
            # println("Average time for $(count_7) trials, time_7 = $(avg_time_7), time_LP_7 = $(avg_LP_7)")
            # string_avg_7 = @sprintf "%.2E" avg_time_7
            # string_avg_LP_7 = @sprintf "%.2E" avg_LP_7
            # string_ppt_LP_7 = @sprintf "%.2E" avg_LP_7 / avg_time_7
            # string_sv_7 = @sprintf "%.2E" avg_sv_7
            # string_cv_7 = @sprintf "%.2E" avg_cv_7
            # string_ppt_sv_7 = @sprintf "%.2f" avg_sv_7 / avg_time_7
            # string_ppt_cv_7 = @sprintf "%.2f" avg_cv_7 / avg_time_7
            # string_ratio_7 = @sprintf "%.2f" avg_time_0 / avg_time_7
            # string_ratio_LP_7 = @sprintf "%.2f" avg_LP_0 / avg_LP_7
            # string_ratio_ppt_LP_7 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_7 / avg_time_7)
        end

        ##################################################
        ##################################################

        if 8 in algos
            info_8 = isf(matrice, options_8)
            count_8 = 1
            
            time_8, FLP_8, ILP_8, time_LP_8, avg_LP_8, prop_LP_8, syms_8, asyms_8, dupli_syms_8, dupli_asyms_8, time_sv_8, prop_sv_8, checks_8, detect_ratio_8, time_cover_8, avg_cover_8, prop_cover_8 = isf_benchmark_values(info_8, options_8)
            
            while count_8 < 3 || (time_8 < 100 && count_8 < 100)
                count_8 += 1
                info_8 = isf(matrice, options_8)
                time_8       += info_8.cput_total
                time_LP_8    += info_8.cput_lp
                prop_LP_8    += (info_8.cput_lp / info_8.cput_total)
                time_sv_8    += info_8.cput_sv
                prop_sv_8    += (info_8.cput_sv / info_8.cput_total)
                time_cover_8 += info_8.cput_cover
                prop_cover_8 += (info_8.cput_cover / info_8.cput_total)
            end
            time_8       /= count_8
            time_LP_8    /= count_8
            prop_LP_8    /= count_8
            time_sv_8    /= count_8
            prop_sv_8    /= count_8
            time_cover_8 /= count_8
            prop_cover_8 /= count_8
            
            push!(DataFrames_vector[DF_index], [time_8, FLP_8, ILP_8, (FLP_0+ILP_0)/(FLP_8+ILP_8) , time_LP_8, time_LP_8 / (FLP_8 + ILP_8), prop_LP_8, syms_8, asyms_8, dupli_syms_8, dupli_asyms_8, time_sv_8, prop_sv_8, checks_8, detect_ratio_8, time_cover_8, 0, prop_cover_8])
            # println("Average time for $(count_8) trials, time_8 = $(avg_time_8), time_LP_8 = $(avg_LP_8)")
            # string_avg_8 = @sprintf "%.2E" avg_time_8
            # string_avg_LP_8 = @sprintf "%.2E" avg_LP_8
            # string_ppt_LP_8 = @sprintf "%.2f" avg_LP_8 / avg_time_8
            # string_ratio_8 = @sprintf "%.2f" avg_time_0 / avg_time_8
            # string_ratio_LP_8 = @sprintf "%.2f" avg_LP_0 / avg_LP_8
            # string_ratio_ppt_LP_8 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_8 / avg_time_8)
        end

        ##################################################
        ##################################################

        if 9 in algos
            info_9 = isf(matrice, options_9)
            count_9 = 1
            
            time_9, FLP_9, ILP_9, time_LP_9, avg_LP_9, prop_LP_9, syms_9, asyms_9, dupli_syms_9, dupli_asyms_9, time_sv_9, prop_sv_9, checks_9, detect_ratio_9, time_cover_9, avg_cover_9, prop_cover_9 = isf_benchmark_values(info_9, options_9)
            
            while count_9 < 3 || (time_9 < 100 && count_9 < 100)
                count_9 += 1
                info_9 = isf(matrice, options_9)
                time_9       += info_9.cput_total
                time_LP_9    += info_9.cput_lp
                prop_LP_9    += (info_9.cput_lp / info_9.cput_total)
                time_sv_9    += info_9.cput_sv
                prop_sv_9    += (info_9.cput_sv / info_9.cput_total)
                time_cover_9 += info_9.cput_cover
                prop_cover_9 += (info_9.cput_cover / info_9.cput_total)
            end
            time_9       /= count_9
            time_LP_9    /= count_9
            prop_LP_9    /= count_9
            time_sv_9    /= count_9
            prop_sv_9    /= count_9
            time_cover_9 /= count_9
            prop_cover_9 /= count_9
            
            push!(DataFrames_vector[DF_index], [time_9, FLP_9, ILP_9, (FLP_0+ILP_0)/(FLP_9+ILP_9) , time_LP_9, time_LP_9 / (FLP_9 + ILP_9), prop_LP_9, syms_9, asyms_9, dupli_syms_9, dupli_asyms_9, time_sv_9, prop_sv_9, checks_9, detect_ratio_9, time_cover_9, 0, prop_cover_9])
            # println("Average time for $(count_9) trials, time_9 = $(avg_time_9), time_LP_9 = $(avg_LP_9)")
            # string_avg_9 = @sprintf "%.2E" avg_time_9
            # string_avg_LP_9 = @sprintf "%.2E" avg_LP_9
            # string_ppt_LP_9 = @sprintf "%.2f" avg_LP_9 / avg_time_9
            # string_ratio_9 = @sprintf "%.2f" avg_time_0 / avg_time_9
            # string_ratio_LP_9 = @sprintf "%.2f" avg_LP_0 / avg_LP_9
            # string_ratio_ppt_LP_9 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_9 / avg_time_9)
        end

        ##################################################
        ##################################################

        if 10 in algos
            info_10 = isf(matrice, options_10)
            count_10 = 1
            
            time_10, FLP_10, ILP_10, time_LP_10, avg_LP_10, prop_LP_10, syms_10, asyms_10, dupli_syms_10, dupli_asyms_10, time_sv_10, prop_sv_10, checks_10, detect_ratio_10, time_cover_10, avg_cover_10, prop_cover_10 = isf_benchmark_values(info_10, options_10)
            
            while count_10 < 3 || (time_10 < 100 && count_10 < 100)
                count_10 += 1
                info_10 = isf(matrice, options_10)
                time_10       += info_10.cput_total
                time_LP_10    += info_10.cput_lp
                prop_LP_10    += (info_10.cput_lp / info_10.cput_total)
                time_sv_10    += info_10.cput_sv
                prop_sv_10    += (info_10.cput_sv / info_10.cput_total)
                time_cover_10 += info_10.cput_cover
                prop_cover_10 += (info_10.cput_cover / info_10.cput_total)
            end
            time_10       /= count_10
            time_LP_10    /= count_10
            prop_LP_10    /= count_10
            time_sv_10    /= count_10
            prop_sv_10    /= count_10
            time_cover_10 /= count_10
            prop_cover_10 /= count_10
            
            push!(DataFrames_vector[DF_index], [time_10, FLP_10, ILP_10, (FLP_0+ILP_0)/(FLP_10+ILP_10) , time_LP_10, time_LP_10 / (FLP_10 + ILP_10), prop_LP_10, syms_10, asyms_10, dupli_syms_10, dupli_asyms_10, time_sv_10, prop_sv_10, checks_10, detect_ratio_10, time_cover_10, 0, prop_cover_10])
            # println("Average time for $(count_10) trials, time_10 = $(avg_time_10), time_LP_10 = $(avg_LP_10)")
            # string_avg_10 = @sprintf "%.2E" avg_time_10
            # string_avg_LP_10 = @sprintf "%.2E" avg_LP_10
            # string_ppt_LP_10 = @sprintf "%.2f" avg_LP_10 / avg_time_10
            # string_ratio_10 = @sprintf "%.2f" avg_time_0 / avg_time_10
            # string_ratio_LP_10 = @sprintf "%.2f" avg_LP_0 / avg_LP_10
            # string_ratio_ppt_LP_10 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_10 / avg_time_10)
        end

        ##################################################
        ##################################################

        if 11 in algos
            info_11 = isf(matrice, options_11)
            count_11 = 1
            
            time_11, FLP_11, ILP_11, time_LP_11, avg_LP_11, prop_LP_11, syms_11, asyms_11, dupli_syms_11, dupli_asyms_11, time_sv_11, prop_sv_11, checks_11, detect_ratio_11, time_cover_11, avg_cover_11, prop_cover_11 = isf_benchmark_values(info_11, options_11)
            
            while count_11 < 3 || (time_11 < 100 && count_11 < 100)
                count_11 += 1
                info_11 = isf(matrice, options_11)
                time_11       += info_11.cput_total
                time_LP_11    += info_11.cput_lp
                prop_LP_11    += (info_11.cput_lp / info_11.cput_total)
                time_sv_11    += info_11.cput_sv
                prop_sv_11    += (info_11.cput_sv / info_11.cput_total)
                time_cover_11 += info_11.cput_cover
                prop_cover_11 += (info_11.cput_cover / info_11.cput_total)
            end
            time_11       /= count_11
            time_LP_11    /= count_11
            prop_LP_11    /= count_11
            time_sv_11    /= count_11
            prop_sv_11    /= count_11
            time_cover_11 /= count_11
            prop_cover_11 /= count_11
            
            push!(DataFrames_vector[DF_index], [time_11, FLP_11, ILP_11, (FLP_0+ILP_0)/(FLP_11+ILP_11) , time_LP_11, time_LP_11 / (FLP_11 + ILP_11), prop_LP_11, syms_11, asyms_11, dupli_syms_11, dupli_asyms_11, time_sv_11, prop_sv_11, checks_11, detect_ratio_11, time_cover_11, 0, prop_cover_11])
            # println("Average time for $(count_11) trials, time_11 = $(avg_time_11), time_LP_11 = $(avg_LP_11)")
            # string_avg_11 = @sprintf "%.2E" avg_time_11
            # string_avg_LP_11 = @sprintf "%.2E" avg_LP_11
            # string_ppt_LP_11 = @sprintf "%.2f" avg_LP_11 / avg_time_11
            # string_ratio_11 = @sprintf "%.2f" avg_time_0 / avg_time_11
            # string_ratio_LP_11 = @sprintf "%.2f" avg_LP_0 / avg_LP_11
            # string_ratio_ppt_LP_11 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_11 / avg_time_11)
        end

        ##################################################
        ##################################################

        if 12 in algos
            info_12 = isf(matrice, options_12)
            count_12 = 1
            
            time_12, FLP_12, ILP_12, time_LP_12, avg_LP_12, prop_LP_12, syms_12, asyms_12, dupli_syms_12, dupli_asyms_12, time_sv_12, prop_sv_12, checks_12, detect_ratio_12, time_cover_12, avg_cover_12, prop_cover_12 = isf_benchmark_values(info_12, options_12)
            
            while count_12 < 3 || (time_12 < 100 && count_12 < 100)
                count_12 += 1
                info_12 = isf(matrice, options_12)
                time_12       += info_12.cput_total
                time_LP_12    += info_12.cput_lp
                prop_LP_12    += (info_12.cput_lp / info_12.cput_total)
                time_sv_12    += info_12.cput_sv
                prop_sv_12    += (info_12.cput_sv / info_12.cput_total)
                time_cover_12 += info_12.cput_cover
                prop_cover_12 += (info_12.cput_cover / info_12.cput_total)
            end
            time_12       /= count_12
            time_LP_12    /= count_12
            prop_LP_12    /= count_12
            time_sv_12    /= count_12
            prop_sv_12    /= count_12
            time_cover_12 /= count_12
            prop_cover_12 /= count_12
            
            push!(DataFrames_vector[DF_index], [time_12, FLP_12, ILP_12, (FLP_0+ILP_0)/(FLP_12+ILP_12) , time_LP_12, time_LP_12 / (FLP_12 + ILP_12), prop_LP_12, syms_12, asyms_12, dupli_syms_12, dupli_asyms_12, time_sv_12, prop_sv_12, checks_12, detect_ratio_12, time_cover_12, time_cover_12 / checks_12, prop_cover_12])
            # println("Average time for $(count_12) trials, time_12 = $(avg_time_12), time_LP_12 = $(avg_LP_12)")
            # string_avg_12 = @sprintf "%.2E" avg_time_12
            # string_avg_LP_12 = @sprintf "%.2E" avg_LP_12
            # string_ppt_LP_12 = @sprintf "%.2E" avg_LP_12 / avg_time_12
            # string_sv_12 = @sprintf "%.2E" avg_sv_12
            # string_cv_12 = @sprintf "%.2E" avg_cv_12
            # string_ppt_sv_12 = @sprintf "%.2f" avg_sv_12 / avg_time_12
            # string_ppt_cv_12 = @sprintf "%.2f" avg_cv_12 / avg_time_12
            # string_ratio_12 = @sprintf "%.2f" avg_time_0 / avg_time_12
            # string_ratio_LP_12 = @sprintf "%.2f" avg_LP_0 / avg_LP_12
            # string_ratio_ppt_LP_12 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_12 / avg_time_12)
            # string_ratio_sv_12 = @sprintf "%.2f" avg_sv_4 / avg_sv_12
            # string_ratio_ppt_sv_12 = @sprintf "%.2f" (avg_sv_4 / avg_time_4) / (avg_sv_12 / avg_time_12)
            # string_ratio_cv_12 = @sprintf "%.2f" avg_cv_4 / avg_cv_12
            # string_ratio_ppt_cv_12 = @sprintf "%.2f" (avg_cv_4 / avg_time_4) / (avg_cv_12 / avg_time_12)
        end

        ##################################################
        ##################################################

        if 13 in algos
            info_13 = isf(matrice, options_13)
            count_13 = 1
            
            time_13, FLP_13, ILP_13, time_LP_13, avg_LP_13, prop_LP_13, syms_13, asyms_13, dupli_syms_13, dupli_asyms_13, time_sv_13, prop_sv_13, checks_13, detect_ratio_13, time_cover_13, avg_cover_13, prop_cover_13 = isf_benchmark_values(info_13, options_13)
            
            while count_13 < 3 || (time_13 < 100 && count_13 < 100)
                count_13 += 1
                info_13 = isf(matrice, options_13)
                time_13       += info_13.cput_total
                time_LP_13    += info_13.cput_lp
                prop_LP_13    += (info_13.cput_lp / info_13.cput_total)
                time_sv_13    += info_13.cput_sv
                prop_sv_13    += (info_13.cput_sv / info_13.cput_total)
                time_cover_13 += info_13.cput_cover
                prop_cover_13 += (info_13.cput_cover / info_13.cput_total)
            end
            time_13       /= count_13
            time_LP_13    /= count_13
            prop_LP_13    /= count_13
            time_sv_13    /= count_13
            prop_sv_13    /= count_13
            time_cover_13 /= count_13
            prop_cover_13 /= count_13
            
            push!(DataFrames_vector[DF_index], [time_13, FLP_13, ILP_13, (FLP_0+ILP_0)/(FLP_13+ILP_13) , time_LP_13, time_LP_13 / (FLP_13 + ILP_13), prop_LP_13, syms_13, asyms_13, dupli_syms_13, dupli_asyms_13, time_sv_13, prop_sv_13, checks_13, detect_ratio_13, time_cover_13, time_cover_13 / checks_13, prop_cover_13])
            # println("Average time for $(count_13) trials, time_13 = $(avg_time_13), time_LP_13 = $(avg_LP_13)")
            # string_avg_13 = @sprintf "%.2E" avg_time_13
            # string_avg_LP_13 = @sprintf "%.2E" avg_LP_13
            # string_ppt_LP_13 = @sprintf "%.2E" avg_LP_13 / avg_time_13
            # string_sv_13 = @sprintf "%.2E" avg_sv_13
            # string_cv_13 = @sprintf "%.2E" avg_cv_13
            # string_ppt_sv_13 = @sprintf "%.2f" avg_sv_13 / avg_time_13
            # string_ppt_cv_13 = @sprintf "%.2f" avg_cv_13 / avg_time_13
            # string_ratio_13 = @sprintf "%.2f" avg_time_0 / avg_time_13
            # string_ratio_LP_13 = @sprintf "%.2f" avg_LP_0 / avg_LP_13
            # string_ratio_ppt_LP_13 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_13 / avg_time_13)
            # string_ratio_sv_13 = @sprintf "%.2f" avg_sv_5 / avg_sv_13
            # string_ratio_ppt_sv_13 = @sprintf "%.2f" (avg_sv_5 / avg_time_5) / (avg_sv_13 / avg_time_13)
            # string_ratio_cv_13 = @sprintf "%.2f" avg_cv_5 / avg_cv_13
            # string_ratio_ppt_cv_13 = @sprintf "%.2f" (avg_cv_5 / avg_time_5) / (avg_cv_13 / avg_time_13)
        end

        ##################################################
        ##################################################

        if 14 in algos
            info_14 = isf(matrice, options_14)
            count_14 = 1
            
            time_14, FLP_14, ILP_14, time_LP_14, avg_LP_14, prop_LP_14, syms_14, asyms_14, dupli_syms_14, dupli_asyms_14, time_sv_14, prop_sv_14, checks_14, detect_ratio_14, time_cover_14, avg_cover_14, prop_cover_14 = isf_benchmark_values(info_14, options_14)
            
            while count_14 < 3 || (time_14 < 100 && count_14 < 100)
                count_14 += 1
                info_14 = isf(matrice, options_14)
                time_14       += info_14.cput_total
                time_LP_14    += info_14.cput_lp
                prop_LP_14    += (info_14.cput_lp / info_14.cput_total)
                time_sv_14    += info_14.cput_sv
                prop_sv_14    += (info_14.cput_sv / info_14.cput_total)
                time_cover_14 += info_14.cput_cover
                prop_cover_14 += (info_14.cput_cover / info_14.cput_total)
            end
            time_14       /= count_14
            time_LP_14    /= count_14
            prop_LP_14    /= count_14
            time_sv_14    /= count_14
            prop_sv_14    /= count_14
            time_cover_14 /= count_14
            prop_cover_14 /= count_14
            
            push!(DataFrames_vector[DF_index], [time_14, FLP_14, ILP_14, (FLP_0+ILP_0)/(FLP_14+ILP_14) , time_LP_14, time_LP_14 / (FLP_14 + ILP_14), prop_LP_14, syms_14, asyms_14, dupli_syms_14, dupli_asyms_14, time_sv_14, prop_sv_14, checks_14, detect_ratio_14, time_cover_14, time_cover_14 / checks_14, prop_cover_14])
            # println("Average time for $(count_14) trials, time_14 = $(avg_time_14), time_LP_14 = $(avg_LP_14)")
            # string_avg_14 = @sprintf "%.2E" avg_time_14
            # string_avg_LP_14 = @sprintf "%.2E" avg_LP_14
            # string_ppt_LP_14 = @sprintf "%.2E" avg_LP_14 / avg_time_14
            # string_sv_14 = @sprintf "%.2E" avg_sv_14
            # string_cv_14 = @sprintf "%.2E" avg_cv_14
            # string_ppt_sv_14 = @sprintf "%.2f" avg_sv_14 / avg_time_14
            # string_ppt_cv_14 = @sprintf "%.2f" avg_cv_14 / avg_time_14
            # string_ratio_14 = @sprintf "%.2f" avg_time_0 / avg_time_14
            # string_ratio_LP_14 = @sprintf "%.2f" avg_LP_0 / avg_LP_14
            # string_ratio_ppt_LP_14 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_14 / avg_time_14)
            # string_ratio_sv_14 = @sprintf "%.2f" avg_sv_6 / avg_sv_14
            # string_ratio_ppt_sv_14 = @sprintf "%.2f" (avg_sv_6 / avg_time_6) / (avg_sv_14 / avg_time_14)
            # string_ratio_cv_14 = @sprintf "%.2f" avg_cv_6 / avg_cv_14
            # string_ratio_ppt_cv_14 = @sprintf "%.2f" (avg_cv_6 / avg_time_6) / (avg_cv_14 / avg_time_14)
        end

        ##################################################
        ##################################################

        if 15 in algos
            info_15 = isf(matrice, options_15)
            count_15 = 1
            
            time_15, FLP_15, ILP_15, time_LP_15, avg_LP_15, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, avg_cover_15, prop_cover_15 = isf_benchmark_values(info_15, options_15)
            
            while count_15 < 3 || (time_15 < 100 && count_15 < 100)
                count_15 += 1
                info_15 = isf(matrice, options_15)
                time_15       += info_15.cput_total
                time_LP_15    += info_15.cput_lp
                prop_LP_15    += (info_15.cput_lp / info_15.cput_total)
                time_sv_15    += info_15.cput_sv
                prop_sv_15    += (info_15.cput_sv / info_15.cput_total)
                time_cover_15 += info_15.cput_cover
                prop_cover_15 += (info_15.cput_cover / info_15.cput_total)
            end
            time_15       /= count_15
            time_LP_15    /= count_15
            prop_LP_15    /= count_15
            time_sv_15    /= count_15
            prop_sv_15    /= count_15
            time_cover_15 /= count_15
            prop_cover_15 /= count_15
            
            push!(DataFrames_vector[DF_index], [time_15, FLP_15, ILP_15, 0, time_LP_15, 0, prop_LP_15, syms_15, asyms_15, dupli_syms_15, dupli_asyms_15, time_sv_15, prop_sv_15, checks_15, detect_ratio_15, time_cover_15, time_cover_15 / checks_15, prop_cover_15])
            # println("Average time for $(count_15) trials, time_15 = $(avg_time_15), time_LP_15 = $(avg_LP_15)")
            # string_avg_15 = @sprintf "%.2E" avg_time_15
            # string_avg_LP_15 = @sprintf "%.2E" avg_LP_15
            # string_ppt_LP_15 = @sprintf "%.2E" avg_LP_15 / avg_time_15
            # string_sv_15 = @sprintf "%.2E" avg_sv_15
            # string_cv_15 = @sprintf "%.2E" avg_cv_15
            # string_ppt_sv_15 = @sprintf "%.2f" avg_sv_15 / avg_time_15
            # string_ppt_cv_15 = @sprintf "%.2f" avg_cv_15 / avg_time_15
            # string_ratio_15 = @sprintf "%.2f" avg_time_0 / avg_time_15
            # string_ratio_LP_15 = @sprintf "%.2f" avg_LP_0 / avg_LP_15
            # string_ratio_ppt_LP_15 = @sprintf "%.2f" (avg_LP_0 / avg_time_0) / (avg_LP_15 / avg_time_15)
            # string_ratio_sv_15 = @sprintf "%.2f" avg_sv_7 / avg_sv_15
            # string_ratio_ppt_sv_15 = @sprintf "%.2f" (avg_sv_7 / avg_time_7) / (avg_sv_15 / avg_time_15)
            # string_ratio_cv_15 = @sprintf "%.2f" avg_cv_7 / avg_cv_15
            # string_ratio_ppt_cv_15 = @sprintf "%.2f" (avg_cv_7 / avg_time_7) / (avg_cv_15 / avg_time_15)
        end

        ##################################################

        if algos == [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15]

            push!(DataFrames_times, [name, time_0, time_8, time_0/time_8, time_1, time_0/time_1, time_9, time_0/time_9, time_2, time_0/time_2, time_10, time_0/time_10, time_3, time_0/time_3, time_11, time_0/time_11, time_4, time_0/time_4, time_12, time_0/time_12, time_5, time_0/time_5, time_13, time_0/time_13, time_6, time_0/time_6, time_14, time_0/time_14, time_7, time_0/time_7, time_15, time_0/time_15])
            ### for the complete set of results
            # println(" & \\multicolumn{4}{c|}{RC} & \\multicolumn{4}{c|}{A} & \\multicolumn{4}{c|}{AB} & \\multicolumn{4}{c|}{ABC}    \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp \\\\ \\hline ")
            # println("time & \\multicolumn{1}{c|}{\$" * string_avg_0 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_8 * "\$} & " * string_ratio_8 * " & \\multicolumn{1}{c|}{\$" * string_avg_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_9 * "\$} & " * string_ratio_9 * " & \\multicolumn{1}{c|}{\$" * string_avg_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_10 * "\$} & " * string_ratio_10 * " & \\multicolumn{1}{c|}{\$" * string_avg_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_11 * "\$} & " * string_ratio_11 * "\\\\ \\hline")
            # println("LP & \\multicolumn{1}{c|}{\$" * string_avg_LP_0 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_8 * "\$} & " * string_ratio_LP_8 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_9 * "\$} & " * string_ratio_LP_9 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_10 * "\$} & " * string_ratio_LP_10 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_11 * "\$} & " * string_ratio_LP_11 * "\\\\ \\hline")
            # println("\\%LP & \\multicolumn{1}{c|}{\$" * string_ppt_LP_0 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_8 * "\$} & " * string_ratio_ppt_LP_8 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_1 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_9 * "\$} & " * string_ratio_ppt_LP_9 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_2 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_10 * "\$} & " * string_ratio_ppt_LP_10 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_3 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_11 * "\$} & " * string_ratio_ppt_LP_11 * "\\\\ \\hline")
            # println("sv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  \\\\ \\hline ")
            # println("\\%sv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  \\\\ \\hline")
            # println("cv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  \\\\ \\hline")
            # println("\\%cv & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} &  & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\multicolumn{1}{c|}{} & \\\\ \\hline ")
            # println("\n")
            # println(" & \\multicolumn{4}{c|}{ABCD1} & \\multicolumn{4}{c|}{ABCD2} & \\multicolumn{4}{c|}{ABCD3}& \\multicolumn{4}{c|}{AD4} \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp  & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp & \\multicolumn{1}{c|}{base}   & \\multicolumn{1}{c|}{comp} & \\multicolumn{1}{c|}{HnH}     & comp \\\\ \\hline ")
            # println("time & \\multicolumn{1}{c|}{\$" * string_avg_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_12 * "\$} & " * string_ratio_12 * " & \\multicolumn{1}{c|}{\$" * string_avg_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_13 * "\$} & " * string_ratio_13 * " & \\multicolumn{1}{c|}{\$" * string_avg_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_14 * "\$} & " * string_ratio_14 * " & \\multicolumn{1}{c|}{\$" * string_avg_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_15 * "\$} & " * string_ratio_15 * "\\\\ \\hline")
            # println("LP & \\multicolumn{1}{c|}{\$" * string_avg_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_12 * "\$} & " * string_ratio_LP_12 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_13 * "\$} & " * string_ratio_LP_13 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_14 * "\$} & " * string_ratio_LP_14 * " & \\multicolumn{1}{c|}{\$" * string_avg_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_avg_LP_15 * "\$} & " * string_ratio_LP_15 * "\\\\ \\hline")
            # println("\\%LP & \\multicolumn{1}{c|}{\$" * string_ppt_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_4 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_12 * "\$} & " * string_ratio_ppt_LP_12 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_5 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_13 * "\$} & " * string_ratio_ppt_LP_13 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_6 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_14 * "\$} & " * string_ratio_ppt_LP_14 * " & \\multicolumn{1}{c|}{\$" * string_ppt_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ratio_ppt_LP_7 * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_LP_15 * "\$} & " * string_ratio_ppt_LP_15 * "\\\\ \\hline")
            # println("sv & \\multicolumn{1}{c|}{\$" * string_sv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_12 * "\$} & " * string_ratio_sv_12 * " & \\multicolumn{1}{c|}{\$" * string_sv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_13 * "\$} & " * string_ratio_sv_13 * "& \\multicolumn{1}{c|}{\$" * string_sv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_14 * "\$} & " * string_ratio_sv_14 * " & \\multicolumn{1}{c|}{\$" * string_sv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_sv_15 * "\$} & " * string_ratio_sv_15 * "\\\\ \\hline")
            # println("\\%sv & \\multicolumn{1}{c|}{\$" * string_ppt_sv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_12 * "\$} & " * string_ratio_ppt_sv_12 * " & \\multicolumn{1}{c|}{\$" * string_ppt_sv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_13 * "\$} & " * string_ratio_ppt_sv_13 * "& \\multicolumn{1}{c|}{\$" * string_ppt_sv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_14 * "\$} & " * string_ratio_ppt_sv_14 * " & \\multicolumn{1}{c|}{\$" * string_ppt_sv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_sv_15 * "\$} & " * string_ratio_ppt_sv_15 * "\\\\ \\hline")
            # println("cv & \\multicolumn{1}{c|}{\$" * string_cv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_12 * "\$} & " * string_ratio_cv_12 * " & \\multicolumn{1}{c|}{\$" * string_cv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_13 * "\$} & " * string_ratio_cv_13 * "& \\multicolumn{1}{c|}{\$" * string_cv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_14 * "\$} & " * string_ratio_cv_14 * " & \\multicolumn{1}{c|}{\$" * string_cv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_cv_15 * "\$} & " * string_ratio_cv_15 * "\\\\ \\hline")
            # println("\\%cv & \\multicolumn{1}{c|}{\$" * string_ppt_cv_4 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_12 * "\$} & " * string_ratio_ppt_cv_12 * " & \\multicolumn{1}{c|}{\$" * string_ppt_cv_5 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_13 * "\$} & " * string_ratio_ppt_cv_13 * "& \\multicolumn{1}{c|}{\$" * string_ppt_cv_6 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_14 * "\$} & " * string_ratio_ppt_cv_14 * " & \\multicolumn{1}{c|}{\$" * string_ppt_cv_7 * "\$} & \\multicolumn{1}{c|}{\$" * " " * "\$} & \\multicolumn{1}{c|}{\$" * string_ppt_cv_15 * "\$} & " * string_ratio_ppt_cv_15 * "\\\\ \\hline")
            ### end of the complete results
            # println("\n")

            # println("Problem & \\multicolumn{2}{c|}{RC} & \\multicolumn{2}{c|}{A} & \\multicolumn{2}{c|}{AB} & \\multicolumn{2}{c|}{ABC} & \\multicolumn{2}{c|}{D1} & \\multicolumn{2}{c|}{D2} & \\multicolumn{2}{c|}{D3} & \\multicolumn{2}{c|}{D4} \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_0 * "\$} & " * " " * " & \\multicolumn{1}{c|}{\$" * string_avg_1 * "\$} & {\\color{blue}{" * string_ratio_1* "}}& \\multicolumn{1}{c|}{\$" * string_avg_2 * "\$} & {\\color{blue}{" * string_ratio_2* "}} & \\multicolumn{1}{c|}{\$" * string_avg_3 * "\$} & {\\color{blue}{" * string_ratio_3* "}}& \\multicolumn{1}{c|}{\$" * string_avg_4 * "\$} & {\\color{blue}{" * string_ratio_4* "}} & \\multicolumn{1}{c|}{\$" * string_avg_5 * "\$} & {\\color{blue}{" * string_ratio_5* "}}& \\multicolumn{1}{c|}{\$" * string_avg_6 * "\$} & {\\color{blue}{" * string_ratio_6* "}} & \\multicolumn{1}{c|}{\$" * string_avg_7 * "\$} & {\\color{blue}{" * string_ratio_7* "}}\\\\ \\hline")
            # println("\n")
            # println("Problem & \\multicolumn{2}{c|}{RC} & \\multicolumn{2}{c|}{A} & \\multicolumn{2}{c|}{AB} & \\multicolumn{2}{c|}{ABC} & \\multicolumn{2}{c|}{D1} & \\multicolumn{2}{c|}{D2} & \\multicolumn{2}{c|}{D3} & \\multicolumn{2}{c|}{D4} \\\\ \\hline")
            # println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_8 * "\$} & {\\color{blue}{" * string_ratio_8 * "}} & \\multicolumn{1}{c|}{\$" * string_avg_9 * "\$} & {\\color{blue}{" * string_ratio_9* "}}& \\multicolumn{1}{c|}{\$" * string_avg_10 * "\$} & {\\color{blue}{" * string_ratio_10* "}} & \\multicolumn{1}{c|}{\$" * string_avg_11 * "\$} & {\\color{blue}{" * string_ratio_11* "}}& \\multicolumn{1}{c|}{\$" * string_avg_12 * "\$} & {\\color{blue}{" * string_ratio_12* "}} & \\multicolumn{1}{c|}{\$" * string_avg_13 * "\$} & {\\color{blue}{" * string_ratio_13* "}}& \\multicolumn{1}{c|}{\$" * string_avg_14 * "\$} & {\\color{blue}{" * string_ratio_14* "}} & \\multicolumn{1}{c|}{\$" * string_avg_15 * "\$} & {\\color{blue}{" * string_ratio_15* "}}\\\\ \\hline")

            # println("\n\n")
        end

    end

    return DataFrames_vector, DataFrames_times

end
   
