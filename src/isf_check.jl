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
`options, info = isf_check(options, values, info)`
preliminary check on some values of options
"""
function isf_check(options::Options, values::Values, info::Info)

    # checking options.algo
    if options.algorithm < 0 || options.algorithm > 15 || !(options.algorithm isa Int)
        info.flag = values.fail_argument
        println("\n### isf_check: options.algorithm should be a nonnegative integer <= 15\n")
    end

    # checking options.bestv
    if options.bestv < 0 || options.bestv > 3 || !(options.bestv isa Int)
        info.flag = values.fail_argument
        println("\n### isf_check: options.bestv should be an integer <= 3\n")
    end

    # checking options.dvnear0
    if !(options.dvnear0 isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.dvnear0 should be an boolean\n")
    end

    # checking options.rational
    if !(options.rational isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.rational should be a boolean\n")
    end

    # checking options.s
    if !(options.s isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.s should be boolean\n")
    end

    # checking options.sv
    if options.sv < 0 || options.sv > 5 || !(options.sv isa Int)
        info.flag = values.fail_argument
        println("\n### isf_check: options.sv should be an integer <= 4\n")
    end

    # checking options.svsprod
    if !(options.svsprod isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.svsprod should be an boolean\n")
    end

    # checking options.symmetry
    if ~(options.symmetry isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.symmetry is not a boolean\n")
    end

    # checking options.wechelon
    if !(options.wechelon isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.wechelon should be an boolean\n")
    end

    # checking options.withd
    if !(options.withd isa Bool)
        info.flag = values.fail_argument
        println("\n### isf_check: options.withd should be boolean\n")
    end

    return options, info
end