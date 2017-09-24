function elem_max(x)
    index = findmax(abs.(x))[2]
    return x[index]
end

function second_largest(x)
    if length(x) == 1
        return x[1]
    end
    return x[2] != x[1] ? x[2] : second_largest(x[2:length(x)])
end
