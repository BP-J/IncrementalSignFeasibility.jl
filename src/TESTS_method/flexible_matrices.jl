using BenchmarkTools
using ElasticArrays
using Random
using LinearAlgebra

"""
Data structure to store vectors in a matrix, adding columns as needed.
When the matrix gets full, adding a new column results in the size doubling.
"""
mutable struct FlexibleMatrixColumnVects
    mat::Matrix{Int64}    # data structure
    ncols::Int64          # number of actually filled columns
end

function FlexibleMatrixColumnVects(nrows)
    FlexibleMatrixColumnVects(zeros(Int64, nrows, 1), 0)
end

function hcat!(A::FlexibleMatrixColumnVects, x::Vector)
    n, p = size(A.mat)
    if p == A.ncols
        A.mat = hcat(A.mat, zeros(Int,n,p))
    end
    A.ncols += 1
    A.mat[:, A.ncols] .= x
    return A
end

function test_flexible_hcat(n,p)
    A = FlexibleMatrixColumnVects(n)
    for _ in 1:p
        hcat!(A, 2 * bitrand(n) .- 1)
    end
    return A
end

function test_normal_hcat(n,p)
    A = Matrix{Int64}(undef, n, 0)
    for _ in 1:p
        A = hcat(A, 2 * bitrand(n) .- 1)
    end
    return A
end

function test_elastic_hcat(n,p)
    A = ElasticMatrix{Int64}(undef, n, 0)
    for _ in 1:p
        append!(A, 2 * bitrand(n) .- 1)
    end
    return A
end


# println("Testing \"flexible\" hcat with 20 rows and 1000 columns")
# @btime test_flexible_hcat(20, 1000);
# println("Testing \"flexible\" hcat with 50 rows and 1000 columns")
# @btime test_flexible_hcat(50, 1000);
# println("Testing \"flexible\" hcat with 50 rows and 2000 columns")
# @btime test_flexible_hcat(50, 2000);
# println("Testing \"flexible\" hcat with 50 rows and 5000 columns")
# @btime test_flexible_hcat(50, 5000);

# println()

# println("Testing \"normal\" hcat with 20 rows and 1000 columns")
# @btime test_normal_hcat(20, 1000);
# println("Testing \"normal\" hcat with 50 rows and 1000 columns")
# @btime test_normal_hcat(50, 1000);
# println("Testing \"normal\" hcat with 50 rows and 2000 columns")
# @btime test_normal_hcat(50, 2000);
# println("Testing \"normal\" hcat with 50 rows and 5000 columns")
# @btime test_normal_hcat(50, 5000);

# println()

# println("Testing \"elastic\" hcat with 20 rows and 1000 columns")
# @btime test_elastic_hcat(20, 1000);
# println("Testing \"elastic\" hcat with 50 rows and 1000 columns")
# @btime test_elastic_hcat(50, 1000);
# println("Testing \"elastic\" hcat with 50 rows and 2000 columns")
# @btime test_elastic_hcat(50, 2000);
# println("Testing \"elastic\" hcat with 50 rows and 5000 columns")
# @btime test_elastic_hcat(50, 5000);

"""
Data structure to store vectors in a matrix, adding rows as needed.
When the matrix gets full, adding a new row results in the size doubling.
"""
mutable struct FlexibleMatrixRowVects
    mat::Matrix{Int64}    # data structure
    nrows::Int64          # number of actually filled rows
end

function FlexibleMatrixRowVects(ncols)
    FlexibleMatrixRowVects(zeros(Float64, 1, ncols), 0)
end

function vcat!(A::FlexibleMatrixRowVects, x::Vector)
    n, p = size(A.mat)
    if n == A.nrows
        A.mat = vcat(A.mat, zeros(Int,n,p))
    end
    A.nrows += 1
    A.mat[A.nrows, :] .= x
    return A
end

function test_flexible_vcat(n,p)
    A = FlexibleMatrixRowVects(p)
    for _ in 1:n
        vcat!(A, 2 * bitrand(p) .- 1)
    end
    return A
end

function test_normal_vcat(n,p)
    A = Matrix{Float64}(undef, 0, p)
    for _ in 1:n
        A = vcat(A, 2 * bitrand(1,p) .- 1)
    end
    return A
end

# function test_elastic_vcat(n,p)
#     A = ElasticMatrix{Float64}(undef, 0, p)
#     for _ in 1:n
#         append!(A, 2 * bitrand(p) .- 1)
#     end
#     return A
# end


# println("Testing \"flexible\" vcat with 20 columns and 100 rows")
# @btime test_flexible_vcat(20, 100);
# println("Testing \"flexible\" vcat with 50 columns and 100 rows")
# @btime test_flexible_vcat(50, 100);
# println("Testing \"flexible\" vcat with 50 columns and 200 rows")
# @btime test_flexible_vcat(50, 200);
# println("Testing \"flexible\" vcat with 50 columns and 500 rows")
# @btime test_flexible_vcat(50, 500);

# println()

# println("Testing \"normal\" vcat with 20 columns and 100 rows")
# @btime test_normal_vcat(20, 100);
# println("Testing \"normal\" vcat with 50 columns and 100 rows")
# @btime test_normal_vcat(50, 100);
# println("Testing \"normal\" vcat with 50 columns and 200 rows")
# @btime test_normal_vcat(50, 200);
# println("Testing \"normal\" vcat with 50 columns and 500 rows")
# @btime test_normal_vcat(50, 500);

# println()

# println("Testing \"elastic\" vcat with 20 columns and 100 rows")
# @btime test_elastic_vcat(20, 100);
# println("Testing \"elastic\" vcat with 50 columns and 100 rows")
# @btime test_elastic_vcat(50, 100);
# println("Testing \"elastic\" vcat with 50 columns and 200 rows")
# @btime test_elastic_vcat(50, 200);
# println("Testing \"elastic\" vcat with 50 columns and 500 rows")
# @btime test_elastic_vcat(50, 500);