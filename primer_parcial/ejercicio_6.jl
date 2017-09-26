include("lib/google_matrix.jl")
include("lib/power_method.jl")

#EJERCICIO 6
n = 10000
A = zeros(2n+1, 2n+1);
q₀ = ones(2n+1)

A[1, :] = A[:, 1] = q₀
for i = 2:2n+1
    if i + 2 <= (2n+1)
        A[i, i+2] = 1
    else
        A[i, (i+2)%(2n+1)+1] = 1
    end
end

Gₐ = crear_matriz_google(A, 0.85);
rₐ = metodo_potencia(Gₐ', q₀, 1e-10, var = "classic")[1]
rₐ = rₐ / norm(rₐ, 1)

#La probabilidad de terminar en la página 0 converge a 0.3.

#Cambiando los aristas por bidireccionales
B = A + A'
B[1, :] = B[:, 1] = q₀
Gᵦ = crear_matriz_google(B, 0.15);
rᵦ = metodo_potencia(Gᵦ', q₀, 1e-10, var = "classic")[1]
rᵦ= rᵦ / norm(rᵦ, 1)

#= La distribución estacionaria no cambia porque intuitivamente, como todos
los estados eran ya positivos recurrentes, no añadimos ninguna posibilidad
de "atorarnos" ni logramos "entrar a algo nuevo". De hecho, la cadena
original ya era irreducible=#
