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
`M = isf_generation_rand_float_affine(dimension, vector_number, bound)`
Returns a matrix in R^{dimension+1 x vector_number}, the last row corresponding 
to the right-hand sides with floating values between -bound and bound.
"""
function isf_generation_rand_float_affine(dimension::Int, vector_number::Int, bound::Float64)

    return bound * (2*rand(dimension+1, vector_number) .- 1.0)
end

"""
`M = isf_generation_rand_int_affine(dimension, vector_number, bound)`
Returns a matrix in R^{dimension+1 x vector_number}, the last row corresponding 
to the right-hand sides with integer values between -bound and bound.
"""
function isf_generation_rand_int_affine(dimension::Int, vector_number::Int, bound::Int)

    return rand(-bound:+bound, dimension+1, vector_number)
end

"""
`M = isf_generation_perm_affine(dimension)`
Returns a matrix corresponding to the affine permutahedron.
"""
function isf_generation_perm_affine(dimension::Int)
    
    M = zeros(Int, dimension+1, Int(dimension*(dimension+1)/2))
    M[dimension+1,1:dimension] = ones(Int, dimension)               # right-hand sides

    indice = dimension
    for i in 1:dimension
        M[i,i] = 1              # identity part 
        for j in i+1:dimension
            indice += 1  # the first dimension ones are taken, so the 
            M[i, indice] = +1
            M[j, indice] = -1
        end
    end
    return M
end

"""
`M = isf_generation_2d_affine(dimension, vector_number, bound)`
Returns a matrix with 2 blocks corresponding to the fan in the upper right 
and the general position in smaller dimension in the bottom left. 
"""
function isf_generation_2d_affine(dimension::Int, vector_number::Int, bound::Int)

    M = zeros(Int, dimension+1, vector_number)
    M[3:dimension+1, 1:dimension-2] = rand(-bound:bound, dimension-1, dimension-2)                          # the lower left part square block
    M[[1,2,dimension+1], dimension-1:vector_number] = rand(-bound:bound, 3, vector_number-dimension+2)    # the upper right 'fan' block
    return M
end

"""
`M = isf_generation_ratio_affine(dimension, vector_numer, bound, ratio)`
Returns a matrix following the definition rules of the 
instances using a ratio: random vectors and rhs form the first dimension ones, 
then each remaining may be random (if rand() > r) or 
a linear combination (weights ∈ [-10, +10]), with the number of vectors 
in the combination as well as the weights being random.
"""
function isf_generation_ratio_affine(dimension::Int, vector_number::Int, bound, ratio::Float64)

    M = zeros(typeof(bound), dimension+1, vector_number)
    M[:,1:dimension] = rand(-bound:bound, dimension+1,dimension)    # the random part

    for i in dimension+1:vector_number
        p = rand()
        if p < ratio
            elements = bitrand(i-1) # element taken if 1, not taken if 0
            while sum(elements) < 2
                elements = bitrand(i-1)     # ensure the combination has at least 2 indices
            end
            println(i)
            println("rand() < ratio")
            println(elements)
            for j in 1:i-1
                if elements[j] == 1
                    coef = rand(-3:3)
                    println(coef)
                    M[:,i] += coef * M[:,j] # linear combination
                else
                    println(0)
                end
            end
        else
            println(i)
            println("rand() > ratio")
            println("random")
            M[:,i] = rand(-bound:bound, dimension+1)    # the random part
        end
    end
    return M
end

"""
`M = isf_generation_srand_affine(dimension, vector_number, bound, number)`
Returns a matrix following the rules of the (affine) srand construction:
the first dimension vectors represent the identity matrix, while 
the remaining vector_number - dimension vectors have number nonzero coordinates
in [-10:+10], and the right-hand sides have values in [-bound:bound]. 
"""
function isf_generation_srand_affine(dimension::Int, vector_number::Int, bound, number)

    M = zeros(typeof(bound), dimension+1, vector_number)
    for i in 1:dimension
        M[i,i] = 1
    end

    coefs = [-10:-1 ; 1:10]

    for i in 1:vector_number-dimension
        vecs_LC = sample(1:dimension, number, replace=false, ordered=true)  # to have no duplicates + in the right order
        for j in 1:number
            M[vecs_LC[j],i+dimension] = rand(coefs)
        end
        if rand() > 0.5 # random right-hand sides
            M[dimension+1,i+dimension] = rand(-bound:bound)
        end
    end
    return M
end


##############################
# How the data was generated
##############################

# rand : ?, ?, 5.0
# perm : just dimension
# degen: ?, 20, 20
# ratio: ?, 20, 50, [0.7;0.9]
# srand: 8, 20, 10, [2;4;6]