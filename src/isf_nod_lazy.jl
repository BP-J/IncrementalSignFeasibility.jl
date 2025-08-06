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
`isf_nod_lazy_H!(p, info)`

Version D5 who does the following operations:
from each stem vector, it creates the sign vectors covering it
then groups all these and does setdiff({± 1}^p, complementary).
"""
function isf_nod_lazy_H!(p::Int64, info::Info)

    ### completion of the stems - since each stem of size l leads to 2^{p-l(+1)} elements of S^c, they are treated individually 
    @suppress begin
        tick()
    end
    for stem in eachcol(info.stems_sym)
        isf_nod_completion_H!(stem, findall(!iszero, stem), p, info) 
    end
    info.cput_cover += tok()

    ### setdiff part: preallocate the big matrix and fill up the lines, alternating 0's and 1'
    # {± 1}^p
    M_base = [1 0]
    M = zeros(Bool, p-1, 2^(p-1))
    for i in 1:p-1
        M[i,:] = repeat(M_base, inner=(1,2^(p-1-i)), outer=(1,2^(i-1)))
    end
    M = [ones(Bool, 1, 2^(p-1)); M]

    # setdiff call
    @suppress begin
        tick()
    end
    info.s = map(v -> 2*v .- 1, setdiff(eachcol(M), info.sclst)) # eachcol(reshape(collect(Iterators.flatten(info.sclst)) ,(p,size(info.sclst,1))))
    info.cput_lp += tok()
    info.ns = size(info.s,1)
end

"""
`isf_nod_lazy_nH!(p, info)`

Version D5 who does the following operations:
from each stem vector, it creates the sign vectors covering it
then groups all these and does setdiff({± 1}^p, complementary).
"""
function isf_nod_lazy_nH!(p::Int64, info::Info)

    ### completion of the stems - since each stem of size l leads to 2^{p-l(+1)} elements of S^c, they are treated individually 
    @suppress begin
        tick()
    end
    for stem in eachcol(info.stems_sym)
        isf_nod_completion_nH_sym!(stem, findall(!iszero, stem), p, info) 
    end

    for stem in eachcol(info.stems_asym)
        isf_nod_completion_nH_asym!(stem, findall(!iszero, stem), p, info) 
    end
    info.cput_cover += tok()

    ### setdiff part: preallocate the big matrix and fill up the lines, alternating 0's and 1'
    # {± 1}^p
    M_base = [1 0]
    M = zeros(Bool, p, 2^p)
    for i in 1:p
        M[i,:] = repeat(M_base, inner=(1,2^(p-i)), outer=(1,2^(i-1)))
    end

    # setdiff call
    @suppress begin
        tick()
    end
    info.s = map(v -> 2*v .- 1, setdiff(eachcol(M), info.sclst)) # eachcol(reshape(collect(Iterators.flatten(info.sclst)) ,(p,size(info.sclst,1))))
    info.cput_lp += tok()
    info.ns = Int(size(info.s,1))
end
