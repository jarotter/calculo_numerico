"""

    elem_max(x)

Regresa el elemento de mayor magnitud de un vector
complejo `x`.

"""
function elem_max(x)
    index = findmax(abs.(x))[2]
    return x[index]
end

"""

    second_largest(x)

Regresa el segundo elemento de mayor magnitud de un vector
complejo `x`. En caso de un empate, regresa el de menor
índice.

"""
function second_largest(x)
    if length(x) == 1
        return x[1]
    end
    return x[2] != x[1] ? x[2] : second_largest(x[2:length(x)])
end
