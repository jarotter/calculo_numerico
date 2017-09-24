include("lib/power_method.jl")
include("lib/aux.jl")

#EJERCICIO 2: CASO ρ = 0
ε = 1e-10

#(a)
A = [1 1 2; -1 9 3; 0 -1 3];
n = size(A)[1]
q₀ = [1; 1; 1];
v̂, λ̂ = metodo_potencia(A, q₀, ε, var = "shift");

#(b)
F = eigfact(A);
ind_min = findmin(F[:values])[2];
λ₁ = F[:values][ind_min];
v₁ = F[:vectors][:, ind_min];
display("Julia calculó $λ₁ y nosotros $λ̂")

#(c)
#=Salvo un fator constante, nuestros vectores propios calculados son iguales,
como podemos observar ||⋅||∞-unitarizando el de Julia=#
display(v₁ / norm(v₁, Inf))
display(v̂)

#Sea m el número de razones que quieren calcularse, por ejemplo,
m = 30
q₁, λ = metodo_potencia(A, q₀, ε, 1, var = "shift");
i = 1
while i <= m
    q₂, λ = metodo_potencia(A, q₁, ε, 1, var = "shift");
    r = norm(q₂ - v̂, Inf) / norm(q₁ - v̂, Inf);
    display("la razón $i es $r")
    q₁ = q₂
    i += 1
end

#(d)
λ₂ = second_largest(sort(F[:values]))
r = abs(λ₁/ λ₂);
display("La razón de convergencia teórica es $r")
#=Efectivamente, las razones en (c) y (d) son iguales salvo error de
aproximación numérica. =#


#EJERCICIO 2: CASO ρ = 3.5
#(a)
v̂, λ̂ = metodo_potencia(A, q₀, ε, ρ = 3.5, var = "shift");

#(b)
F = eigfact(A - 3.5*eye(n));
ind_min = findmin(abs.(F[:values]))[2];
λ₁ = F[:values][ind_min] + 3.5;
v₁ = F[:vectors][:, ind_min];
display("Julia calculó $λ₁ y nosotros $λ̂")

#(c)
#=Salvo un fator constante, nuestros vectores propios calculados son iguales,
como podemos observar ||⋅||∞-unitarizando el de Julia=#
display(v₁ / norm(v₁, Inf))
display(v̂)

#Sea m el número de razones que quieren calcularse, por ejemplo,
m = 30
q₁, λ = metodo_potencia(A, q₀, ε, 1, ρ = 3.5, var = "shift");
i = 1
while i <= m
    q₂, λ = metodo_potencia(A, q₁, ε, 1, ρ = 3.5, var = "shift");
    r = norm(q₂ - v̂, Inf) / norm(q₁ - v̂, Inf);
    display("La razón $i es $r")
    q₁ = q₂
    i += 1
end

#(d)
λ₂ = second_largest(sort(abs.(F[:values])))
r = abs((λ₁-3.5)/ λ₂);
display("La razón de convergencia teórica es $r")
#=Efectivamente, las razones en (c) y (d) son iguales salvo error de
aproximación numérica. =#
