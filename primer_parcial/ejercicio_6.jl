#Pkg.add("Distributions")
include("lib/google_matrix.jl")
include("lib/power_method.jl")

#EJERCICIO 6
n = 100
A = zeros(2n+1, 2n+1);
q₀ = ones(2n+1)

A[1, :] = A[:, 1] = q₀
A[1,1]=0
for i = 2:2n+1
    if i + 2 <= (2n+1)
        A[i, i+2] = 1
    else
        A[i, (i+2)%(2n+1)+1] = 1
    end
end

Gₐ = crear_matriz_google(A, 0.15);
rₐ = metodo_potencia(Gₐ', q₀, 1e-10, var = "classic")[1]
rₐ = rₐ / norm(rₐ, 1)
#=Damos pagerank en norma uno para representar probabilidades en la distribución
estacionaria de la cadena de Markov representada por Gₐ=#

#La probabilidad de terminar en la página 0 converge a 0.3.

#Cambiando los aristas por bidireccionales
B = A + A'
B[1, :] = B[:, 1] = q₀
B[1,1]=0
Gᵦ = crear_matriz_google(B, 0.15);
rᵦ = metodo_potencia(Gᵦ', q₀, 1e-10, var = "classic")[1]
rᵦ= rᵦ / norm(rᵦ, 1)
