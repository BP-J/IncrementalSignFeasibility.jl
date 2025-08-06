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
# Not a particular function but defines some basic pieces reused later
# such as structures, type converting functions and so on.

"""
`A = echelonned_form!(A)`

Computes the echelonned form without reduction. 
99.9% copied from the RowEchelon package, but with increased tolerance. 
"""
function echelonned_form!(A::Matrix{T}, ɛ=T <: Union{Rational,Integer} ? 0 : eps(norm(A,Inf))) where T
    nr, nc = size(A)
    i = j = 1
    while i <= nr && j <= nc
        (m, mi) = findmax(abs.(A[i:nr,j]))
        mi = mi+i - 1
        if m <= 100000 * ɛ
            if ɛ > 0
                A[i:nr,j] .= zero(T)
            end
            j += 1
        else
            for k=j:nc
                A[i, k], A[mi, k] = A[mi, k], A[i, k]
            end
            d = A[i,j]
            for k = j:nc
                A[i,k] /= d
            end
            for k = 1:nr
                if k != i
                    d = A[k,j]
                    for l = j:nc
                        A[k,l] -= d*A[i,l]
                    end
                end
            end
            i += 1
            j += 1
        end
    end
    A
end

echelonned_formconv(::Type{T}, A::Matrix) where {T} = echelonned_form!(copyto!(similar(A, T), A))

echelonned_form(A::Matrix{T}) where {T} = echelonned_form!(copy(A))
echelonned_form(A::Matrix{T}) where {T <: Complex} = echelonned_formconv(ComplexF64, A)
echelonned_form(A::Matrix{ComplexF64}) = echelonned_form!(copy(A))
echelonned_form(A::Matrix{T}) where {T <: Union{Integer, Float16, Float32}} = echelonned_formconv(Float64, A)
echelonned_form(A::AbstractMatrix) = echelonned_form(Matrix(A))