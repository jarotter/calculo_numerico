using(Distributions)

function crear_matriz_adyacencia(n, p)
    A = zeros((n,n))
    b = Bernoulli(p)e
    for i = 1:n
        A[i, 1:n] = rand(b, n)
        A[i, i] = 0
    end
    return A
end

function crear_matriz_google(A, q)
    n = size(A)[1]
    K = A * ones(n) #Las sumas por renglón

    for j = 1:n
        A[j, 1:n] = K[j] == 0 ? ones(n)/n : A[j, 1:n] / K[j]
    end

    S = ones((n,n))/n

    G = q*S + (1-q)*A
    return G
end
