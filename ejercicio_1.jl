include("lib/power_method.jl")
include("lib/aux.jl")


#EJERCICIO 1
ε = 1e-10

#(a)
A = [1 1 2; -1 9 3; 0 -1 3];
q₀ = [1; 1; 1];
v̂, λ̂ = metodo_potencia(A, q₀, ε, var = "classic");

#(b)
F = eigfact(A);
ind_max = findmax(F[:values])[2];
λ₁ = F[:values][ind_max];
v₁ = F[:vectors][:, ind_max];
display("Julia calculó $λ₁ y nosotros $λ̂")

#(c)
#=Salvo un fator constante, nuestros vectores propios calculados son iguales,
como podemos observar ||⋅||∞-unitarizando el de Julia=#
display(v₁ / norm(v₁, Inf))
display(v̂)

#Sea m el número de razones que quieren calcularse, por ejemplo,
m = 30
q₁, λ = metodo_potencia(A, q₀, ε, 1, var = "classic");
i = 1
while i <= m
    q₂, λ = metodo_potencia(A, q₁, ε, 1, var = "classic");
    r = norm(q₂ - v̂, Inf) / norm(q₁ - v̂, Inf);
    display("la razón $i es $r")
    q₁ = q₂
    i += 1
end

#(d)
λ₂ = second_largest(sort(F[:values], rev = true));
r = abs(λ₂ / λ₁);
display("La razón de convergencia teórica es $r")
#=Efectivamente, las razones en (c) y (d) son iguales salvo error de
aproximación numérica. =#
