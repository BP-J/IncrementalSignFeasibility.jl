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
`isf_coherence!(options)`
Checks the coherence of the options, puts flag to 1 if not. 
"""
function isf_coherence!(options::Options, values::Values)

    # symmetry Checks
    if options.symmetry && options.algorithm >= 8
        println("Symmetric case with asymetric algorithm")
        info.flag = values.fail_argument
        return
    end

    if options.dvnear0 && (options.algorithm % 8 <= 1)
        println("Modification B requires algorithms 2-7 or 10-15")
        info.flag = values.fail_argument
        return
    end

    if options.bestv > 0 && (options.algorithm % 8 <= 2)
        println("Modification C requires algorithms 3-7 or 11-15")
        info.flag = values.fail_argument
        return
    end

    if !options.withd && (options.algorithm % 8 != 7)
        println("Working without directions requires algorithms 7 or 15")
        info.flag = values.fail_argument
        return
    end

    if options.sv > 0 && (options.algorithm % 8 <= 3)
        println("Versions with stem vectors require algorithms 4-7 or 12-15")
        info.flag = values.fail_argument
        return
    end

    if options.sv == 0 && options.svsprod
        println("Cannot have svsprod true without stem vectors")
        info.flag = values.fail_argument
        return
    end
end