#####
# This file is used to compare the use of 'unique' in the stem vectors structure(s). 

# One option is to store the (potential) stem vectors then use unique to remove duplicates. 
# The other is to check at every addition if the stem vector is already known. 

# Some functions call 'unique' while others manually check ; it is very slow, 
# unique should always be used.

# The second part only uses "unique" after the computation is done. 

#####

include("comparisons_sv_structure.jl")

function comparison_matrix_test(N, p)
    rng = MersenneTwister(0)
    M = zeros(Int, N, p)
    indices = []
    for i in 1:N
        v = 2 * bitrand(rng,p) .- 1
        if !(v in eachrow(M[1:i,:]))
            M[i,:] = v
            push!(indices,i)
        end
    end
    return M[indices,:]
end

# ----------------------------------------

function comparison_matrix_unique(N, p)
    rng = MersenneTwister(0)
    M = zeros(Int, N, p)
    for i in 1:N
        M[i,:] = 2 * bitrand(rng, p) .- 1
    end
    return unique(M,dims=1)
end

# ----------------------------------------

function comparison_list_test(N, p)
    rng = MersenneTwister(0)
    L = []
    for i in 1:N
        v = 2 * bitrand(rng,p) .- 1
        if !(v in L)
            push!(L, v)
        end
    end
    return L
end

# ----------------------------------------

function comparison_list_unique(N, p)
    rng = MersenneTwister(0)
    L = []
    for i in 1:N
        push!(L, 2 * bitrand(rng,p) .- 1)
    end
    return unique(L)
end

# ----------------------------------------

function comparison_inc_matrix_test(N, p)
    rng = MersenneTwister(0)
    M = zeros(Int, 100, p)
    c = 0
    for i in 1:N
        v = 2 * bitrand(rng,p) .- 1
        if !(v in eachrow(M))
            c += 1
            M[c,:] = v
            if c == size(M,1)
                M = [M ; zeros(Int, size(M,1), p)]
            end
        end
    end
    return M[1:c,:]
end

# ----------------------------------------

function comparison_inc_matrix_unique(N, p)
    rng = MersenneTwister(0)
    M = zeros(Int, 100, p)
    for i in 1:N
        M[i,:] = 2 * bitrand(rng,p) .- 1
        if i == size(M,1)
            M = [M ; zeros(Int, size(M,1), p)]
        end
    end
    return unique(M[1:N,:], dims=1)
end

# ----------------------------------------

function comparison_flex_matrix_test(N, p)
    rng = MersenneTwister(0)
    M = FlexibleMatrixRowVects(p)
    c = 0
    for _ in 1:N
        v = 2 * bitrand(rng,p) .- 1
        if !(v in eachrow(M.mat))
            vcat!(M, v)
            c += 1
        end
    end
    return M
end

# ----------------------------------------

function comparison_flex_matrix_unique(N, p)
    rng = MersenneTwister(0)
    M = FlexibleMatrixRowVects(p)
    for _ in 1:N
        vcat!(M, 2 * bitrand(rng, p) .- 1)
    end
    MM = unique(M.mat, dims=1)
    return FlexibleMatrixRowVects(MM, size(MM,1))
end

# ----------------------------------------

function comparison_elastic_test(N, p)
    rng = MersenneTwister(0)
    E = ElasticMatrix{Int64}(undef, p, 0) # transposed later

    for _ in 1:N
        v = 2 * bitrand(rng,p) .- 1
        if !(v in eachcol(E))
            append!(E, v)
        end
    end
    return E
end

# ----------------------------------------

function comparison_elastic_unique(N, p)
    rng = MersenneTwister(0)
    E = ElasticMatrix{Int64}(undef, p, 0) # transposed later

    for _ in 1:N
        append!(E, 2 * bitrand(rng, p) .- 1)
    end
    return unique(E, dims=2)
end

# ----------------------------------------
# ----------------------------------------
# ----------------------------------------

function unique_matrix(M::Matrix)
    return unique(M, dims=1)
end

function unique_list(L::Vector{Any})
    return unique(L)
end

function unique_converted_list(LC::Adjoint)
    return unique(LC, dims=1)
end

function unique_flex(F::FlexibleMatrixRowVects)
    MM = unique(F.mat, dims=1)
    return FlexibleMatrixRowVects(MM, size(MM,1))
end

function unique_elastic(E::ElasticMatrix)
    return unique(E, dims=2)
end



# ----------------------------------------



### cannot even be finished for the 'test' ones
# N = 50000
# p = 20
# @btime comparison_matrix_test($N,$p)
# @btime comparison_matrix_unique($N,$p)
# @btime comparison_list_test($N,$p)
# @btime comparison_list_unique($N,$p)
# @btime comparison_inc_matrix_test($N,$p)
# @btime comparison_inc_matrix_unique($N,$p)
# @btime comparison_flex_matrix_test($N,$p)
# @btime comparison_flex_matrix_unique($N,$p)
# @btime comparison_elastic_test($N,$p)
# @btime comparison_elastic_unique($N,$p)

# N = 100000
# p = 20
# @benchmark comparison_matrix_test($N,$p)
# @benchmark comparison_matrix_unique($N,$p)
# @benchmark comparison_list_test($N,$p)
# @benchmark comparison_list_unique($N,$p)
# @benchmark comparison_inc_matrix_test($N,$p)
# @benchmark comparison_inc_matrix_unique($N,$p)
# @benchmark comparison_flex_matrix_test($N,$p)
# @benchmark comparison_flex_matrix_unique($N,$p)
# @benchmark comparison_elastic_test($N,$p)
# @benchmark comparison_elastic_unique($N,$p)