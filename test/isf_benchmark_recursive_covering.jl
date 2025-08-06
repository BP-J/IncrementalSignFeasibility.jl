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

using DelimitedFiles
using Printf

# for name in list_instances

#     fullname = "../../complementarity-internship-2021/papers/1-bdiffmin-affine/matlab/isf/data/data_" * name
#     matrice = readdlm(fullname)
    
#     println(name)

#     info_D1     = isf(matrice, options_8_4)
#     time_D1 = info_D1.cput_total
#     time_cover_D1 = info_D1.cput_cover
#     list_times_D1 = [time_D1]
#     list_cover_D1 = [time_cover_D1]
#     count_D1 = 1
#     while time_D1 < 100 && count_D1 < 100       # maximum 100 trials
#         count_D1 += 1
#         info_D1 = isf(matrice, options_8_4)
#         time_D1 += info_D1.cput_total
#         time_cover_D1 += info_D1.cput_cover
#         append!(list_times_D1, info_D1.cput_total)
#         append!(list_cover_D1, info_D1.cput_cover)
#     end
#     avg_time_D1 = time_D1 / count_D1
#     avg_time_cover_D1 = time_cover_D1 / count_D1
#     println("Average time for $(count_D1) trials, time_D1 = $(avg_time_D1), time_cover_D1 = $(avg_time_cover_D1)")
#     println(list_times_D1)
#     println(list_cover_D1)
#     string_avg_D1 = @sprintf "%.2E" avg_time_D1
#     string_avg_cover_D1 = @sprintf "%.2E" avg_time_cover_D1

#     ##################################################

#     info_D1_rec = isf(matrice, options_8_4r)
#     time_D1_rec = info_D1_rec.cput_total
#     time_cover_D1_rec = info_D1_rec.cput_cover
#     list_times_D1_rec = [time_D1_rec]
#     list_cover_D1_rec = [time_cover_D1_rec]
#     count_D1_rec = 1
#     while time_D1_rec < 100 && count_D1_rec < 100       # maximum 100 trials
#         count_D1_rec += 1
#         info_D1_rec = isf(matrice, options_8_4r)
#         time_D1_rec += info_D1_rec.cput_total
#         time_cover_D1_rec += info_D1_rec.cput_cover
#         append!(list_times_D1_rec, info_D1_rec.cput_total)
#         append!(list_cover_D1_rec, info_D1_rec.cput_cover)
#     end
#     avg_time_D1_rec = time_D1_rec / count_D1_rec
#     avg_time_cover_D1_rec = time_cover_D1_rec / count_D1_rec
#     println("Average time for $(count_D1_rec) trials, time_D1_rec = $(avg_time_D1_rec), time_cover_D1_rec = $(avg_time_cover_D1_rec)")
#     println(list_times_D1_rec)
#     println(list_cover_D1_rec)
#     println("ratio of avg total time = $(avg_time_D1 / avg_time_D1_rec)")
#     println("ratio of avg cover time = $(avg_time_cover_D1 / avg_time_cover_D1_rec)")
#     string_avg_D1_rec = @sprintf "%.2E" avg_time_D1_rec
#     string_avg_cover_D1_rec = @sprintf "%.2E" avg_time_cover_D1_rec
#     string_ratio_D1 = @sprintf "%.2f" avg_time_D1 / avg_time_D1_rec
#     string_ratio_cover_D1 = @sprintf "%.2f" avg_time_cover_D1 / avg_time_cover_D1_rec

#     ##################################################
#     ##################################################

#     info_D2     = isf(matrice, options_8_5)
#     time_D2 = info_D2.cput_total
#     time_cover_D2 = info_D2.cput_cover
#     list_times_D2 = [time_D2]
#     list_cover_D2 = [time_cover_D2]
#     count_D2 = 1
#     while time_D2 < 100 && count_D2 < 100       # maximum 100 trials
#         count_D2 += 1
#         info_D2 = isf(matrice, options_8_5)
#         time_D2 += info_D2.cput_total
#         time_cover_D2 += info_D2.cput_cover
#         append!(list_times_D2, info_D2.cput_total)
#         append!(list_cover_D2, info_D2.cput_cover)
#     end
#     avg_time_D2 = time_D2 / count_D2
#     avg_time_cover_D2 = time_cover_D2 / count_D2
#     println("Average time for $(count_D2) trials, time_D2 = $(avg_time_D2), time_cover_D2 = $(avg_time_cover_D2)")
#     println(list_times_D2)
#     println(list_cover_D2)
#     string_avg_D2 = @sprintf "%.2E" avg_time_D2
#     string_avg_cover_D2 = @sprintf "%.2E" avg_time_cover_D2

#     ##################################################

#     info_D2_rec = isf(matrice, options_8_5r)

#     time_D2_rec = info_D2_rec.cput_total
#     time_cover_D2_rec = info_D2_rec.cput_cover
#     list_times_D2_rec = [time_D2_rec]
#     list_cover_D2_rec = [time_cover_D2_rec]
#     count_D2_rec = 1
#     while time_D2_rec < 100 && count_D2_rec < 100       # maximum 100 trials
#         count_D2_rec += 1
#         info_D2_rec = isf(matrice, options_8_5r)
#         time_D2_rec += info_D2_rec.cput_total
#         time_cover_D2_rec += info_D2_rec.cput_cover
#         append!(list_times_D2_rec, info_D2_rec.cput_total)
#         append!(list_cover_D2_rec, info_D2_rec.cput_cover)
#     end
#     avg_time_D2_rec = time_D2_rec / count_D2_rec
#     avg_time_cover_D2_rec = time_cover_D2_rec / count_D2_rec
#     println("Average time for $(count_D2_rec) trials, time_D2_rec = $(avg_time_D2_rec), time_cover_D2_rec = $(avg_time_cover_D2_rec)")
#     println(list_times_D2_rec)
#     println(list_cover_D2_rec)
#     println("ratio of avg total time = $(avg_time_D2 / avg_time_D2_rec)")
#     println("ratio of avg cover time = $(avg_time_cover_D2 / avg_time_cover_D2_rec)")
#     string_avg_D2_rec = @sprintf "%.2E" avg_time_D2_rec
#     string_avg_cover_D2_rec = @sprintf "%.2E" avg_time_cover_D2_rec
#     string_ratio_D2 = @sprintf "%.2f" avg_time_D2 / avg_time_D2_rec
#     string_ratio_cover_D2 = @sprintf "%.2f" avg_time_cover_D2 / avg_time_cover_D2_rec

#     ##################################################
#     ##################################################

#     info_D3     = isf(matrice, options_8_6w)
#     time_D3 = info_D3.cput_total
#     time_cover_D3 = info_D3.cput_cover
#     list_times_D3 = [time_D3]
#     list_cover_D3 = [time_cover_D3]
#     count_D3 = 1
#     while time_D3 < 100 && count_D3 < 100       # maximum 100 trials
#         count_D3 += 1
#         info_D3 = isf(matrice, options_8_6w)
#         time_D3 += info_D3.cput_total
#         time_cover_D3 += info_D3.cput_cover
#         append!(list_times_D3, info_D3.cput_total)
#         append!(list_cover_D3, info_D3.cput_cover)
#     end
#     avg_time_D3 = time_D3 / count_D3
#     avg_time_cover_D3 = time_cover_D3 / count_D3
#     println("Average time for $(count_D3) trials, time_D3 = $(avg_time_D3), time_cover_D3 = $(avg_time_cover_D3)")
#     println(list_times_D3)
#     println(list_cover_D3)
#     string_avg_D3 = @sprintf "%.2E" avg_time_D3
#     string_avg_cover_D3 = @sprintf "%.2E" avg_time_cover_D3

#     ##################################################

#     info_D3_rec = isf(matrice, options_8_6wr)

#     time_D3_rec = info_D3_rec.cput_total
#     time_cover_D3_rec = info_D3_rec.cput_cover
#     list_times_D3_rec = [time_D3_rec]
#     list_cover_D3_rec = [time_cover_D3_rec]
#     count_D3_rec = 1
#     while time_D3_rec < 100 && count_D3_rec < 100       # maximum 100 trials
#         count_D3_rec += 1
#         info_D3_rec = isf(matrice, options_8_6wr)
#         time_D3_rec += info_D3_rec.cput_total
#         time_cover_D3_rec += info_D3_rec.cput_cover
#         append!(list_times_D3_rec, info_D3_rec.cput_total)
#         append!(list_cover_D3_rec, info_D3_rec.cput_cover)
#     end
#     avg_time_D3_rec = time_D3_rec / count_D3_rec
#     avg_time_cover_D3_rec = time_cover_D3_rec / count_D3_rec
#     println("Average time for $(count_D3_rec) trials, time_D3_rec = $(avg_time_D3_rec), time_cover_D3_rec = $(avg_time_cover_D3_rec)")
#     println(list_times_D3_rec)
#     println(list_cover_D3_rec)
#     println("ratio of avg total time = $(avg_time_D3 / avg_time_D3_rec)")
#     println("ratio of avg cover time = $(avg_time_cover_D3 / avg_time_cover_D3_rec)")
#     string_avg_D3_rec = @sprintf "%.2E" avg_time_D3_rec
#     string_avg_cover_D3_rec = @sprintf "%.2E" avg_time_cover_D3_rec
#     string_ratio_D3 = @sprintf "%.2f" avg_time_D3 / avg_time_D3_rec
#     string_ratio_cover_D3 = @sprintf "%.2f" avg_time_cover_D3 / avg_time_cover_D3_rec

#     ##################################################
#     ##################################################

#     info_D4     = isf(matrice, options_8_7w)
#     time_D4 = info_D4.cput_total
#     time_cover_D4 = info_D4.cput_cover
#     list_times_D4 = [time_D4]
#     list_cover_D4 = [time_cover_D4]
#     count_D4 = 1
#     while time_D4 < 100 && count_D4 < 100       # maximum 100 trials
#         count_D4 += 1
#         info_D4 = isf(matrice, options_8_7w)
#         time_D4 += info_D4.cput_total
#         time_cover_D4 += info_D4.cput_cover
#         append!(list_times_D4, info_D4.cput_total)
#         append!(list_cover_D4, info_D4.cput_cover)
#     end
#     avg_time_D4 = time_D4 / count_D4
#     avg_time_cover_D4 = time_cover_D4 / count_D4
#     println("Average time for $(count_D4) trials, time_D4 = $(avg_time_D4), time_cover_D4 = $(avg_time_cover_D4)")
#     println(list_times_D4)
#     println(list_cover_D4)
#     string_avg_D4 = @sprintf "%.2E" avg_time_D4
#     string_avg_cover_D4 = @sprintf "%.2E" avg_time_cover_D4

#     ##################################################

#     info_D4_rec = isf(matrice, options_8_7wr)

#     time_D4_rec = info_D4_rec.cput_total
#     time_cover_D4_rec = info_D4_rec.cput_cover
#     list_times_D4_rec = [time_D4_rec]
#     list_cover_D4_rec = [time_cover_D4_rec]
#     count_D4_rec = 1
#     while time_D4_rec < 100 && count_D4_rec < 100       # maximum 100 trials
#         count_D4_rec += 1
#         info_D4_rec = isf(matrice, options_8_7wr)
#         time_D4_rec += info_D4_rec.cput_total
#         time_cover_D4_rec += info_D4_rec.cput_cover
#         append!(list_times_D4_rec, info_D4_rec.cput_total)
#         append!(list_cover_D4_rec, info_D4_rec.cput_cover)
#     end
#     avg_time_D4_rec = time_D4_rec / count_D4_rec
#     avg_time_cover_D4_rec = time_cover_D4_rec / count_D4_rec
#     println("Average time for $(count_D4_rec) trials, time_D4_rec = $(avg_time_D4_rec), time_cover_D4_rec = $(avg_time_cover_D4_rec)")
#     println(list_times_D4_rec)
#     println(list_cover_D4_rec)
#     println("ratio of avg total time = $(avg_time_D4 / avg_time_D4_rec)")
#     println("ratio of avg cover time = $(avg_time_cover_D4 / avg_time_cover_D4_rec)")
#     string_avg_D4_rec = @sprintf "%.2E" avg_time_D4_rec
#     string_avg_cover_D4_rec = @sprintf "%.2E" avg_time_cover_D4_rec
#     string_ratio_D4 = @sprintf "%.2f" avg_time_D4 / avg_time_D4_rec
#     string_ratio_cover_D4 = @sprintf "%.2f" avg_time_cover_D4 / avg_time_cover_D4_rec

#     ##################################################
#     ##################################################
#     println()
#     println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_D1 * "\$} & \$" * string_avg_cover_D1 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_D1_rec * "\$} & \$" * string_avg_cover_D1_rec * "\$ & \\multicolumn{1}{c|}{" * string_ratio_D1 * "} & " * string_ratio_cover_D1 * "\\\\ \\hline")
#     println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_D2 * "\$} & \$" * string_avg_cover_D2 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_D2_rec * "\$} & \$" * string_avg_cover_D2_rec * "\$ & \\multicolumn{1}{c|}{" * string_ratio_D2 * "} & " * string_ratio_cover_D2 * "\\\\ \\hline")
#     println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_D3 * "\$} & \$" * string_avg_cover_D3 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_D3_rec * "\$} & \$" * string_avg_cover_D3_rec * "\$ & \\multicolumn{1}{c|}{" * string_ratio_D3 * "} & " * string_ratio_cover_D3 * "\\\\ \\hline")
#     println("\\texttt{$(name)} & \\multicolumn{1}{c|}{\$" * string_avg_D4 * "\$} & \$" * string_avg_cover_D4 * "\$ & \\multicolumn{1}{c|}{\$" * string_avg_D4_rec * "\$} & \$" * string_avg_cover_D4_rec * "\$ & \\multicolumn{1}{c|}{" * string_ratio_D4 * "} & " * string_ratio_cover_D4 * "\\\\ \\hline")

# end