module ExactLinearAlgebra

function echelon(M::Matrix{N}) where N <: Number
    M_aug = vcat(M, eye(N, size(M,2)))
    i,j = size(M)
    while i > 0 && j > 0
        k = findprev(M_aug[i,:],j)
        if k == 0
            i -= 1
            continue
        elseif k < j
            M_aug[:,[k,j]] = M_aug[:,[j,k]]
        end

        @assert !iszero(M_aug[i,j])

        for k=(j-1):-1:1
            if !iszero(M_aug[i,k])
                M_aug[:,k] = M_aug[i,j] * M_aug[:,k] - M_aug[i,k] * M_aug[:,j]
            end
        end
        i -= 1
        j -= 1
    end

    return M_aug

end

function colspan(M::AbstractMatrix{N}) where N <: Number
    M_aug = echelon(M)

    nonzero_cols = [ j for j in indices(M, 2) if any(m != 0 for m in M_aug[indices(M,1),j])]

    return M_aug[indices(M,1), nonzero_cols]
end

function kernel(M::AbstractMatrix{N}) where N <: Number
    M_aug = echelon(M)

    zero_cols = [ j for j in indices(M, 2) if all(m == 0 for m in M_aug[indices(M,1),j])]

    return M_aug[(size(M,1)+1):end, zero_cols]
end


end
