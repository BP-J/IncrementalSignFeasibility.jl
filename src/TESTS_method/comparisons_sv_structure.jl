#####
# These few files are used to know which structure seems better, mainly for the stem vectors. 

# - one option is to store them in a list of vectors, then checking one by one with the current stem vector. 
# - another option is to create a matrix then increase its size every time
# - - or to increase some of the time 

# The Random package is used for the following: 
# - fixing the rng seed (necessary?)
# - using the bitrand function, to obtain 'sign' vectors (then 2 * v .- 1)

# the "launching code" is at the end


## For instance, setting a rng seed and using it makes so that bitrand will be consistent. 

#####

using Random
using LinearAlgebra
using BenchmarkTools

include("flexible_matrices.jl")

function test_construction_liste(N, p)
    L = Any[]
    rng = MersenneTwister(0)

    for i in 1:N
        # v = 2 * bitrand(rng, p) .- 1 # only -1's and +1's
        push!(L, 2 * bitrand(rng, p) .- 1)
    end
    # done row by row as it is the way the stem vectors are added

    return L
end

# ----------------------------------------

function test_construction_matrice(N,p)
    rng = MersenneTwister(0)
    M = zeros(Int, N, p)
    for i in 1:N
        M[i,:] = 2 * bitrand(rng,p) .- 1
    end
    # done row by row as it is the way the stem vectors are added
    return M
end

# ----------------------------------------

function test_construction_increasing_matrice(N, p)

    rng = MersenneTwister(0)
    
    # constructed then row by row
    M = zeros(Int, 100, p)

    for i in 1:N
        M[i,:] = 2 * bitrand(rng,p) .- 1
        if i == size(M,1)
            M = [M ; zeros(Int, size(M,1), p)]
        end
    end

    return M[1:N,:]
end

# ----------------------------------------

function test_construction_flexible_matrix(N, p)

    rng = MersenneTwister(0)

    M = FlexibleMatrixRowVects(p)
    for _ in 1:N
        vcat!(M, 2 * bitrand(rng, p) .- 1)
    end
    return M
end

# ----------------------------------------

function test_construction_elastic_arrays(N, p)
    
    rng = MersenneTwister(0)

    M = ElasticMatrix{Int64}(undef, p, 0) # transposed later

    for _ in 1:N
        append!(M, 2 * bitrand(rng, p) .- 1)
    end
    return M
end

# ----------------------------------------

function conversion_liste(L)
    return hcat(L...)'      # does the right thing
end


# ----------------------------------------
# ----------------------------------------
# ----------------------------------------

# construction of the elements
# N = ...
# p = ...

# L = test_construction_liste(N, p);
# M = test_construction_matrice(N, p);
# F = test_construction_flexible_matrix(N, p);
# E = test_construction_elastic_arrays(N, p);

# @benchmark test_construction_liste($N,$p)
# @benchmark test_construction_matrice($N, $p)
# @benchmark test_construction_increasing_matrice($N, $p)
# @benchmark test_construction_flexible_matrix($N, $p)
# @benchmark test_construction_elastic_arrays($N, $p)
# @benchmark conversion_liste($L)

