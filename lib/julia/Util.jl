import Base: findnext
function findnext{S <: SparseVector}(v::S, i::Int)
    n = searchsortedfirst(v.nzind, i)
    if n > length(v.nzind)
        return 0
    else
        return v.nzind[n]
    end
end

function findnext(m::SparseMatrixCSC, i::Int)
   if i > length(m)
       return 0
   end
   row, col = ind2sub(m, i)
   lo, hi = m.colptr[col], m.colptr[col+1]
   n = searchsortedfirst(m.rowval[lo:hi-1], row)
   if n == 0
       return sub2ind(m, m.rowval[lo], col)
   end
   if n <= hi-lo
       return sub2ind(m, m.rowval[lo+n-1], col)
   end
   nextcol = findnext(c->(c>hi), m.colptr, col+1)
   if nextcol == 0
       return 0
   end
   nextlo = m.colptr[nextcol-1]
   return sub2ind(m, m.rowval[nextlo], nextcol-1)
end


