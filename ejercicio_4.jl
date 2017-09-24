#Ejercicio 4

A = [1 -4 -6; -12 -8 -6; 11 10 10]
F = eigfact(A);
abs.(F[:values])

#=Así pues, no funciona porque no hay un valor propio dominante; el de mayor
magnitud es complejo y por tanto viene con su conjugado, que tiene la misma
magnitud =#
