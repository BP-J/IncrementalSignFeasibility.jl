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
# Computes the upper bound for the central and noncentral
# cases from the number of vectors p and the rank r.
# Both correspond to the "general position" configuration, though
# the definition of general position depends on the case.

"""
`w = isf_central_max(p, r)`
The central formula, attained in (linear) general position. 
The number of vectors is p and the rank r.
"""
function isf_central_max(p::Int64, r::Int64)
    pm = p-1
    smax = 0
    for i in 0:r-1
        smax += binomial(pm, i)
    end
    return 2*smax
end

"""
`z = isf_noncentral_max(p, r)`
The noncentral formula, attained in (affine) general position. 
The number of vectors is p and the rank r. 
"""
function isf_noncentral_max(p::Int64, r::Int64)
    smax = 0
    for i in 0:r
        smax += binomial(p, i)
    end
    return smax
end
