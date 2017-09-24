function qr_simple(A, k, ε)
    A = hessfact(A)[:H]
    while diag(A, -1) < ε
        Q, R = qr(A)
        A = R*Q
    end
    return diag(A)
end

function qr_dynamic(A, k, ε)
    n = size(A)[1]
    λ = zeros(n)
    k = n
    A = hessfact(A)[:H]
    while k > 1
        while abs(A[k, k-1]) > ε
            ρ = A[k, k]
            Q, R = qr(A - ρ*eye(n))
            A = R*Q + ρ*eye(n)
        end
        λ[k] = A[k, k]
        k -= 1
        A = A[1:k, 1:k]
    end
    return λ
end
