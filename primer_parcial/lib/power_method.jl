include("aux.jl")

"""

    metodo_potencia(A, q₀, ε, max_iter = 100; var = "rayleigh", ρ = 0)

Función genérica variaciones del método de la potencia. En todas, el
parámetro `A` es la matriz sobre la cuál se trabaja,
desde un vector inicial `q₀`, con tolerancia ε para
calcular su valor propio principal.

Pueden controlarse también las iteraciones con
`max_iter`.

El argumento `var` puede recibir "rayleigh" para el método
de la potencia inversa con shift de Rayleigh, "classic"
para el tradicional o "shift" para la potencia inversa
con shift (el shift dado en el argumento `ρ`).
"""
function metodo_potencia(A, q₀, ε, max_iter = 100; var = "rayleigh", ρ = 0)
    if var == "classic"
        _potencia_classic(A, q₀, max_iter, ε)
    elseif var == "shift"
        _potencia_inversa_con_shift(A, q₀, max_iter, ε, ρ)
    else
        _potencia_inversa_rayleigh(A, q₀, max_iter)
    end
end

function _potencia_classic(A, q₀, max_iter, ε)
    q = q̂ = q₀
    σ = 0.0

    for j = 1:max_iter
        q = A * q̂
        σ = elem_max(q)
        q = q / σ
        if norm(q - q̂, Inf) < ε
            break
        end
         q̂ = q
    end

    return q, σ
end

function _potencia_inversa_con_shift(A, q₀, max_iter, ε, ρ)
    q = q̂ = q₀
    σ = 0.0
    n = size(A)[1]
    A = A - ρ*eye(n)

    for j = 1:max_iter
        q = A\q̂
        σ = elem_max(q)
        q = q / σ
        if norm(q - q̂, Inf) < ε
            break
        end
        q̂ = q
    end

    return q, ((1 / σ) + ρ)
end

function _potencia_inversa_rayleigh(A, q₀, max_iter)
    q = q̂ = q₀
    n = size(A)[1];
    σ = ρ = 0
    for i = 1:max_iter
        ρ = (q̂' * A * q̂)/(q̂'*q̂)
        Â = A - ρ * eye(n)
        if rank(Â) < n
            break
        end
        q̂ = Â \ q̂
        σ = elem_max(q̂)
        q̂ = q̂ / σ
        q = q̂
    end

    return q, ρ
end
