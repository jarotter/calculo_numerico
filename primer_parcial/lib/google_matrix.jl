using(Distributions)

"""

    crear_matriz_adyacencia(n, p)

Crea una matriz de adyacencia de `n` × `n` donde cada
entrada es uno o cero según un experimento Bernoulli de
parámetro `p`.

"""
function crear_matriz_adyacencia(n, p)
    A = zeros((n,n))
    b = Bernoulli(p)e
    for i = 1:n
        A[i, 1:n] = rand(b, n)
        A[i, i] = 0
    end
    return A
end

"""

    crear_matriz_google(A, q)

A partir de una matriz de adyacencia `A`, crea la matriz
de Google necesaria para ejecutar el algoritmo Pagerank.
La matriz de Google que usamos es irreducible,
aperiódica y estocástica por filas, y
representa las probabilidades de pasar de la página i a la
j pensado como estados en una cadena de Markov homogénea.

El parámetro `q` compensa por el evento probable de que
en cualquier momento el usuario cambie a una página aleatoria
y no sólo a una con hipervínculo desde la actual.
"""

function crear_matriz_google(A, q)
    n = size(A)[1]
    B = A
    K = B * ones(n) #Las sumas por renglón

    for j = 1:n
        B[j, 1:n] = K[j] == 0 ? ones(n)/n : B[j, 1:n] / K[j]
    end

    S = ones((n,n))/n

    G = q*S + (1-q)*B
    return G
end
