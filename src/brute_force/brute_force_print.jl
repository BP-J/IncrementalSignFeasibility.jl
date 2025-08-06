### 
# flag = brute_force_print(V, bvec, d, val, ns, options, values)

# print bvec + verification if needed. Does not depend on the version.

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

function brute_force_print(bvec, d, val, ns, options, success)
    flag = success

    if options.verb2 >= 2
        if options.fout2 == 1
            println("$(ns) | ")
        else
            open(options.fout2, "a") do f
                write(f, "\n$(ns) | ")
            end
        end
    end

    if options.fout2 == 1
        println("$(arr_to_str_bis(bvec))")
    else
        open(options.fout2, "a") do f
            write(f, "\n$(arr_to_str_bis(bvec))")
        end
    end

    if options.verb2 >= 2
        if options.fout2 == 1
            println(" | $(val)")
            println(" |")
            println(" $(d)")
        else
            open(options.fout2, "a") do f
                write(f, " | $(val)")
                write(f, " |")
                write(f, " $(d)")
            end
        end
    end

    return flag
end