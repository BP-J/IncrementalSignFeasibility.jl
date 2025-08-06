### 
# ns, info = brute_force_nH(V, options)

# Brute force algorithm listing the sign vectors by checking every 
# possible system, using a linear optimization problem. See section 5.2.1 (?)
# of the paper mentioned below. 
# Asymmetric version, for matrices \tilde{V} = [V;T] with T = (τ_i)_i. 
# The structure is very similar

include("../isf_tools.jl")
include("brute_force_minus.jl")
include("brute_force_print.jl")
include("brute_force_optim.jl")

#-----------------------------------------------------------------------
#-----------------------------------------------------------------------
# Version 0.1, May, 2023.
#
# If you found this piece of software useful, please cite the paper:
# J.-P. Dussault, J.Ch. Gilbert, B. Plaquevent-Jourdain,
# '?????', 2023.
#
# Authors:
# - Jean-Pierre Dussault (Univ. of Sherbrooke, Canada),
# - Jean Charles Gilbert (INRIA, France),
# - Baptiste Plaquevent-Jourdain (INRIA & Univ. of Sherbrooke, Canada).
#
# Copyright 2023, INRIA (France) and Univ. of Sherbrooke (Canada).
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
# <http://doc.trolltech.com/3.0/license.html>.
#
#-----------------------------------------------------------------------
#-----------------------------------------------------------------------

function brute_force_nH(Vt, options)
    
    # time shenanigans
    date = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")

    success = Int(0)
    ns = 0
    info = isf_get_info()
    info.flag = success

    eps2 = 100*eps()        # tolerance

    n = size(Vt,1) - 1
    p = size(Vt,2)

    if options.verb >= 2
        dline = "--------------------------------------------------"
        eline = "=================================================="
        tline = "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        open(options.fout, "a") do f
            write(f, eline)
        end
    end

    if options.verb >= 1
        if options.fout == 1
            print("\nISF (Version 0.1, 2023-05)")
            print("\n^^^^^^^^^^^^^^^^^^^^^^^^^^")
            print("\nDate / time: ", date)
            print("\nVector dimension (n) = $(n)")
            print("\nNumber of vectors (p) = $(p)")
        else
            open(options.fout, "a") do f
                write(f, "\nISF (Version 0.1, 2023-05)")
                write(f, "\n^^^^^^^^^^^^^^^^^^^^^^^^^^")
                write(f, "\nDate / time: ", date)
                write(f, "\nVector dimension (n) = $(n)")
                write(f, "\nNumber of vectors (p) = $(p)")
            end
        end   # end of options.fout == 1 / == "name..."
    end     # end of options.verb == 1

    # brute force method
    
    # loop on the 2^p nodes: bvec = (s+1)/2 is the binary representation of s = (2*bvec)-1

    bvec = ones(Int, p)  # initial binary vector of ones, by convention
    ns = 0

    # main loop on half the binary vectors: by symmetry one can consider simply the ones for which 
    # bveec[1] = 1

    while true # ?
        
        # main sign vector
        VS= (2*bvec .- 1)' .* Vt
        d, val = brute_force_optim_nH(VS)
        if val < 0
            ns += 1
            push!(info.s, copy(bvec))
            
            if options.verb2 == 1
                if options.fout2 == 1
                    println(arr_to_str_bis(bvec))
                else
                    open(options.fout2, "a") do f
                        write(f, "\n$(arr_to_str_bis(bvec))")
                    end
                end
            end
            if options.verb2 == 2
                info.flag = brute_force_print(bvec, d, val, ns, options, success)
                if info.flag == 1
                    return ns, info
                end
            end
        end

        if options.verb2 >= 3
            info.flag = brute_force_print(bvec, d, val, ns, options, success)
            if info.flag == 1
                return ns, info
            end
        end

        ### part with the opposite sign vector

        
        op_bvec = opposite(bvec)
        VS= (2*op_bvec .- 1)' .* Vt
        d, val = brute_force_optim_nH(VS)
        if val < 0
            ns += 1
            push!(info.s, copy(op_bvec))
            
            if options.verb2 == 1
                if options.fout2 == 1
                    println(arr_to_str_bis(op_bvec))
                else
                    open(options.fout2, "a") do f
                        write(f, "\n$(arr_to_str_bis(op_bvec))")
                    end
                end
            end
            if options.verb2 == 2
                info.flag = brute_force_print(op_bvec, d, val, ns, options, success)
                if info.flag == 1
                    return ns, info
                end
            end
        end

        if options.verb2 >= 3
            info.flag = brute_force_print(op_bvec, d, val, ns, options, success)
            if info.flag == 1
                return ns, info
            end
        end

        # next binary vector

        b = brute_force_minus(bvec[2:p])

        if isempty(b)
            break
        end
        bvec[2:p] = b

    end

    if options.verb >= 1
        if options.fout == 1
            println("Number of sign vectors = $(ns)")
            # time shenanigans
        else
            open(options.fout, "a") do f
                write(f, "\nNumber of sign vectors = $(ns)")
            end
        end
    end

    return ns, info

end