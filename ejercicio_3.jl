include("lib/power_method.jl")

#EJERCICIO 3
ε = 1e-10
A = [1 1 2; -1 9 3; 0 -1 3];
n = size(A)[1]
Λ = ones(n) * NaN
V = zeros((n,n))
guessed = 0
q₀ = ones(n)

while guessed < n
    q₀ = rand(n)
    v, λ = metodo_potencia(A, q₀, ε)
    if norm(Λ - λ, Inf) > 10e-8 #Es linealmente independiente
        guessed += 1
        Λ[guessed] = λ
        V[:, guessed] = v
    else #Hay que revisar independencia lineal
        real_vecs = hcat(V[1:n, 1:guessed], v)'
        if rank(real_vecs) == guessed + 1
            guessed += 1
            Λ[guessed] = λ
            V[:, guessed] = v
        end
    end
end

display(V)
display(Λ)

#=Para verificar que la convergencia es cuadrática, replicamos lo  hecho
en el ejercicio 1, de nuevo siendo m el número de razones a calcular =#
m = 30
v = V[1:n, n]
q₁, λ = metodo_potencia(A, q₀, ε);
i = 1
while i <= m
    q₂, λ = metodo_potencia(A, q₁, ε);
    r = norm(q₂ - v, Inf) / norm(q₁ - v, Inf)^2;
    display("la razón $i es $r")
    q₁ = q₂
    i += 1
end
