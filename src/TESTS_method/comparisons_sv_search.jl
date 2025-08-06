



function build_tests_in(M, n)
    """
    M is the given matrix of 'data'
    size(M) = N, (p)
    n is the number of test elements that are from M

    The function takes some (vaguely equally spaced) lines of the matrix
    """
    N = size(M,1)
    step = Int(floor(N/(n+1)))

    indices = zeros(Int, n)
    for i in 1:n
        indices[i] = i*step + rand(-step:step)
    end
    return M[indices,:]
end

# ----------------------------------------

function build_tests_out(M, n)
    """
    M is the given matrix of 'data'
    size(M) = N, (p)
    n is the number of test elements that are from M

    the lines of the matrix are then made false by adding a zero to ensure they will not be recognized
    """
    N = size(M,1)
    step = Int(floor(N/(n+1)))

    indices = zeros(Int, n)
    for i in 1:n
        indices[i] = i*step + rand(-step:step)
    end
    MM = M[indices,:]
    MM[:,1] = zeros(Int, n) # ±1 changed to 0 to ensure the vectors will not be recognized
    return MM
end

# ----------------------------------------

function multiple_research_matrix(M, MM, n, p)

    for i in 1:n
        V = M * MM[i,:]
        stop = findfirst(x -> abs(x) == p, V)
    end
end

# ----------------------------------------

function multiple_research_list(L, MM, n, p)

    for i in 1:n
        for j in eachindex(L)
            if abs(dot(L[j], MM[i,:])) == p
                break
            end
        end
        # V = [dot(L[j], MM[i,:]) for j in eachindex(L)]
        # stop = findfirst(x -> abs(x) == p, V)
    end
end

# ----------------------------------------

function multiple_research_elastic(E, MM, n, p)

    for i in 1:n
        V = E' * MM[i,:]
        stop = findfirst(x -> abs(x) == p, V)
    end
end

# ----------------------------------------

function recherche_unique_for_liste(L, p, v)
    stop = 0
    # dummy test
    for i in 1:eachindex(L)
        if abs(dot(L[i], v)) == p
            stop = i
            break
        end
    end

    return stop
end

# ----------------------------------------

function recherche_unique_all_liste(L, p, v)
    stop = 0
    
    # complete construction
    Values = [abs(dot(L[i], v)) for i in 1:eachindex(L)]
    stop = findfirst(x -> x == p, Values)
    return stop
end

# ----------------------------------------

function recherche_unique_conversion(L, p, v)
    M = conversion_liste(L)

    stop = 0
    V = M * v
    stop = findfirst(x -> abs(x) == p, V)

    return stop
end

# ----------------------------------------

function recherche_unique_matrice(M, p, v)

    stop = 0
    # complete construction
    V = M * v
    stop = findfirst(x -> abs(x) == p, V)

    return stop
end

# ----------------------------------------

function recherche_unique_flexible(M, p, v)

    V = M.mat * v
    stop = findfirst(x -> abs(x) == p, V)
    return stop
end

# ----------------------------------------

function recherche_unique_elastic(M, p, v)

    V = M' * v
    stop = findfirst(x -> abs(x) == p, V)
    return stop
end


# @benchmark recherche_unique_for_liste($L, $p, $v)
# @benchmark recherche_unique_all_liste($LISTE, $p, $v)
# @benchmark recherche_unique_conversion($LISTE, $p, $v)
# @benchmark recherche_unique_matrice($M, $p, $v)
# @benchmark recherche_unique_flexible($F, $p, $v)
# @benchmark recherche_unique_elastic($E, $p, $v)