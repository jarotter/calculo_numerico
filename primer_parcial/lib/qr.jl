""""

    qr_simple(A, max_iter, ε)

Método QR simple para aproximar valores propios de `A`, con un error
máximo de `ε`, en a lo más `max_iter` iteraciones.

Regresa:

*Un  vector`Λ`, cuyas entradas son los valores propios de `A`.

*Un valor lógico para indicar si `A` es simétrica o no.

*Una matriz `V` con los vectores propios de `A` (sólo cuando `A` es simétrica.
Si no lo es, puede ignorarse)

*El número de iteraciones del algoritmo antes de la convergencia.


"""
function qr_simple(A, max_iter, ε)

    is_symmetric = issymmetric(A)
    F = hessfact(A)
    A = F[:H]
    V = F[:Q]
    num_iter = 0

    for i = 1:max_iter
        Q, R = qr(A)
        V = V*Q
        A = R*Q
        if norm(diag(A, -1)) < ε
            num_iter = i
            break
        end
    end
    return diag(A), is_symmetric, V, min(num_iter, max_iter)
end


""""

    qr_dynamic(A, num_vaps, ε)

Método QR con shift de Wilkinson para aproximar valores propios de `A`,
con un error máximo de `ε`.  Usar el shift de Wilkinson asegura que
el algoritmo siempre converge, cuadráticamente en el peor de los casos, e
incluso puede con valores propios simétricos (a diferencia del shift
con cociente de Rayleigh).
Si no se quieren todos los valores propios,
pueden calcularse sólo `num_vaps`.

Regresa:

*Un vector `Λ`, cuyas entradas son los valoresvpropios de `A`.

*El número de iteraciones del algoritmo antes de la convergencia.


"""
function qr_dynamic(A, num_vaps, ε)

    k = min(num_vaps, size(A)[1])
    λ = zeros(k)
    F = hessfact(A)
    A = F[:H]
    num_iter = 0

    while k > 1

        while abs(A[k, k-1]) > ε
            #Calculamos por partes el shift de Wilkinson
            n = size(A)[1]
            B = A[n-1:n, n-1:n]
            δ = (B[1, 1] - B[2, 2])/2
            δ̂ = δ == 0 ? 1 : sign(δ)
            ρ = B[2, 2] - δ̂*B[2, 1]^2/(abs(δ) + sqrt(δ^2+ B[2,1]^2))

            Q, R = qr(A - ρ*eye(size(A)[1]))
            A = R*Q + ρ*eye(size(A)[1])
            num_iter += 1
        end

        λ[k] = A[k, k]
        k -= 1
        A = A[1:k, 1:k]
    end
    λ[1] = A[1, 1]

    return λ, num_iter
end
