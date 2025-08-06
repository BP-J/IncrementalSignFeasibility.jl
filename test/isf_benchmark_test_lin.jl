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

# types of instances 
list_instances_rand = ["rand_4_8_2.txt", "rand_7_8_4.txt", "rand_7_9_4.txt", "rand_7_10_5.txt", "rand_7_11_4.txt", "rand_7_12_6.txt", "rand_7_13_5.txt", "rand_7_14_7.txt", "rand_8_15_7.txt", "rand_9_16_8.txt", "rand_10_17_9.txt"]
list_instances_2d = ["degen2d_20_4.txt", "degen2d_20_5.txt", "degen2d_20_6.txt", "degen2d_20_7.txt", "degen2d_20_8.txt"]
list_instances_srand = ["srand_8_20_2.txt", "srand_8_20_4.txt", "srand_8_20_6.txt"]
list_instances_perm = ["perm_5.txt", "perm_6.txt", "perm_7.txt", "perm_8.txt"]
list_instances_ratio = ["ratio_20_3_7.txt", "ratio_20_3_9.txt", "ratio_20_4_7.txt", "ratio_20_4_9.txt", "ratio_20_5_7.txt", "ratio_20_5_9.txt", "ratio_20_6_7.txt", "ratio_20_6_9.txt", "ratio_20_7_7.txt", "ratio_20_7_9.txt"]

list_instances_threshold = ["threshold_3.txt", "threshold_4.txt", "threshold_5.txt"]
list_instances_resonance = ["resonance_4.txt", "resonance_5.txt"] # , "resonance_6.txt" # resonance_6 is excluded "a priori" since the computation of its stem vectors takes too long. 
list_instances_demicube = ["demicube_4.txt", "demicube_5.txt", "demicube_6.txt"]
list_instances_crosspoly = ["crosspoly_5.txt", "crosspoly_6.txt", "crosspoly_7.txt", "crosspoly_8.txt", "crosspoly_9.txt", "crosspoly_10.txt", "crosspoly_11.txt", "crosspoly_12.txt"]

# example with all the instances except resonance_6
list_instances = [list_instances_rand ; list_instances_2d ; list_instances_srand ; list_instances_perm ; list_instances_ratio ; list_instances_threshold ; list_instances_resonance ; list_instances_demicube ; list_instances_crosspoly]

# explicit possible options
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

# example with all the tree algorithms
algos=["0","1","2","3","4","4r","5","5r","6","6r","6wr","6w","7","7r","7wr","7w"] # ,"l" is excluded a priori since it may be very slow

# Depending on how many algorithms and how many instances, the rest of the code may be very long. 

for name in list_instances

    fullname = "linear_data/" * name
    matrice = readdlm(fullname)
    println(name, "\n")

    # This part computes the different information for each algorithm on a given instance. 
    if "0" in algos
        info_0 = isf(matrice, options_0)
        println("fin 0,", info_0.cput_total)
    end
    
    if "1" in algos
        info_1 = isf(matrice, options_1)
        println("fin 1,", info_1.cput_total)
    end
    
    if "2" in algos
        info_2 = isf(matrice, options_2)
        println("fin 2,", info_2.cput_total)
    end
    
    if "3" in algos
        info_3 = isf(matrice, options_3)
        println("fin 3,", info_3.cput_total)
    end
    
    if "4" in algos
        info_4 = isf(matrice, options_4)
        println("fin 4,", info_4.cput_total)
    end
    
    if "4r" in algos
        info_4r = isf(matrice, options_4r)
        println("fin 4r,", info_4r.cput_total)
    end
    
    if "5" in algos
        info_5 = isf(matrice, options_5)
        println("fin 5,", info_5.cput_total)
    end
    
    if "5r" in algos
        info_5r = isf(matrice, options_5r)
        println("fin 5r,", info_5r.cput_total)
    end
    
    if "6" in algos
        info_6 = isf(matrice, options_6)
        println("fin 6,", info_6.cput_total)
    end
    
    if "6r" in algos
        info_6r = isf(matrice, options_6r)
        println("fin 6r,", info_6r.cput_total)
    end
    
    if "6wr" in algos
        info_6wr = isf(matrice, options_6wr)
        println("fin 6wr,", info_6wr.cput_total)
    end
    
    if "6w" in algos
        info_6w = isf(matrice, options_6w)
        println("fin 6w,", info_6w.cput_total)
    end
    
    if "7" in algos
        info_7 = isf(matrice, options_7)
        println("fin 7,", info_7.cput_total)
    end
    
    if "7r" in algos
        info_7r = isf(matrice, options_7r)
        println("fin 7r,", info_7r.cput_total)
    end
    
    if "7wr" in algos
        info_7wr = isf(matrice, options_7wr)
        println("fin 7wr,", info_7wr.cput_total)
    end
    
    if "7w" in algos
        info_7w = isf(matrice, options_7w)
        println("fin 7w,", info_7w.cput_total)
    end
    
    if "l" in algos
        info_l = isf(matrice, options_l)
        println("fin l,", info_l.cput_total)
    end
    
    ##################################################
    ##################################################
    ##################################################
    ##################################################

    # quick verification of the number of chambers
    println(info_0.ns, " ", info_1.ns, " ", info_2.ns, " ", info_3.ns, " ", info_4.ns, " ", info_4r.ns, " ", info_5.ns, " ", info_5r.ns, " ", info_6.ns, " ", info_6r.ns, " ", info_6wr.ns, " ", info_6w.ns, " ", info_7.ns, " ", info_7r.ns, " ", info_7wr.ns, " ", info_7w.ns) # , " ", info_l.ns

    # comparisons ; the order is to try comparing the most similar sets (which should accelerate the comparison)
    # to adjust if not all algorithms are to be compared
    
    if info_0.ns != info_1.ns
        # println(length(setdiff(info_0.s, info_1.s)))
        # println(length(setdiff(info_1.s, info_0.s)))
        println("$(name) a un problème entre 0 et 1")
    else                                                # this one is different because the non-RC versions have modification A, which changes the order
        sd = setdiff(info_0.s, [info_1.s ; -info_1.s])  # whereas the RC version, in the symmetric case, has the first sign always +1
        if length(sd) == 0
            println("il y a bien égalité entre 0 et 1")
        else
            println(length(sd))
            println("$(name) a un problème entre 0 et 1")
        end
    end

    if info_1.ns != info_2.ns
        # println(length(setdiff(info_1.s, info_2.s)))
        # println(length(setdiff(info_2.s, info_1.s)))
        println("$(name) a un problème entre 1 et 2")
    else
        sd = setdiff(info_1.s, info_2.s)
        if length(sd) == 0
            println("il y a bien égalité entre 1 et 2")
        else
            println(length(sd))
            println("$(name) a un problème entre 1 et 2")
        end
    end

    if info_2.ns != info_3.ns
        # println(length(setdiff(info_2.s, info_3.s)))
        # println(length(setdiff(info_3.s, info_2.s)))
        println("$(name) a un problème entre 2 et 3")
    else
        sd = setdiff(info_2.s, info_3.s)
        if length(sd) == 0
            println("il y a bien égalité entre 2 et 3")
        else
            println(length(sd))
            println("$(name) a un problème entre 2 et 3")
        end
    end

    if info_3.ns != info_4.ns
        # println(length(setdiff(info_3.s, info_4.s)))
        # println(length(setdiff(info_4.s, info_3.s)))
        println("$(name) a un problème entre 3 et 4")
    else
        sd = setdiff(info_3.s, info_4.s)
        if length(sd) == 0
            println("il y a bien égalité entre 3 et 4")
        else
            println(length(sd))
            println("$(name) a un problème entre 3 et 4")
        end
    end

    if info_4.s != info_4r.s
        println("$(name) a un problème aux algorithmes 4?")
    else
        println("il y a bien égalité entre 4 et 4r")
    end

    if info_4.ns != info_5.ns
        # println(length(setdiff(info_4.s, info_5.s)))
        # println(length(setdiff(info_5.s, info_4.s)))
        println("$(name) a un problème entre 4 et 5")
    else
        sd = setdiff(info_4.s, info_5.s)
        if length(sd) == 0
            println("il y a bien égalité entre 4 et 5")
        else
            println(length(sd))
            println("$(name) a un problème entre 4 et 5")
        end
    end

    if info_5.s != info_5r.s
        println("$(name) a un problème aux algorithmes 5?")
    else
        println("il y a bien égalité entre 5 et 5r")
    end

    if info_5.ns != info_6.ns
        # println(length(setdiff(info_5.s, info_6.s)))
        # println(length(setdiff(info_6.s, info_5.s)))
        println("$(name) a un problème entre 5 et 6")
    else
        sd = setdiff(info_5.s, info_6.s)
        if length(sd) == 0
            println("il y a bien égalité entre 5 et 6")
        else
            println(length(sd))
            println("$(name) a un problème entre 5 et 6")
        end
    end

    if (info_6.s != info_6r.s) || (info_6.s != info_6wr.s) || (info_6.s != info_6w.s) 
        println("$(name) a un problème aux algorithmes 6?")
    else
        println("il y a bien égalité entre 6, 6r, 6rw, 7w")
    end

    if info_6.ns != info_7.ns
        # println(length(setdiff(info_6.s, info_7.s)))
        # println(length(setdiff(info_7.s, info_6.s)))
        println("$(name) a un problème entre 6 et 7")
    else
        sd = setdiff(info_6.s, info_7.s)
        if length(sd) == 0
            println("il y a bien égalité entre 6 et 7")
        else
            println(length(sd))
            println("$(name) a un problème entre 6 et 7")
        end
    end

    if (info_7.s != info_7r.s) || (info_7.s != info_7wr.s) || (info_7.s != info_7w.s) 
        println("$(name) a un problème aux algorithmes 7?")
    else
        println("il y a bien égalité entre 7, 7r, 7rw, 7w")
    end

    if info_7.ns != info_l.ns
        # println(length(setdiff(info_7.s, info_l.s)))
        # println(length(setdiff(info_l.s, info_7.s)))
        println("$(name) a un problème entre 7 et l")
    else                                                # this one is different because the non-RC versions have modification A, which changes the order
        sd = setdiff(info_7.s, [info_l.s ; -info_l.s])                # whereas this version, in the symmetric case, has no specific order
        if length(sd) == 0
            println("il y a bien égalité entre 7 et l")
        else
            println(length(sd))
            println("$(name) a un problème entre 7 et l")
        end
    end

    println("Done")
end
