module Util

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

type BoundedPriorityQueue{K,V}
    values::Vector{Pair{K,V}}
    cur_length::Int
    max_length::Int
    BoundedPriorityQueue(max_length::Int) = new(resize!(Pair{K,V}[], max_length+1), 0, max_length)
end

using Base.Collections: heapparent, percolate_down!, percolate_up!
function bounded_heapify!(xs::AbstractArray, limit, o::Base.Ordering)
    for i in heapparent(limit):-1:1
        percolate_down!(xs, i, o, limit)
    end
    xs
end

function enqueue!{K,V}(x::BoundedPriorityQueue{K,V}, k::K, v::V)
    x.values[x.cur_length+1] = (k=>v)
    if(x.cur_length < x.max_length)
        x.cur_length +=1
    end
    percolate_up!(x.values, x.cur_length, Base.Order.By(a->a[2]))
end

function dequeue!{K,V}(x::BoundedPriorityQueue{K,V})
    if x.cur_length < 1
        throw(BoundsError())
    end
    k,v = x.values[1]
    x.values[1] = x.values[x.cur_length]
    x.cur_length -= 1
    percolate_down!(x.values, 1, Base.Order.By(a->a[2]), x.cur_length)
    return k
end

function peek{K,V}(x::BoundedPriorityQueue{K,V})
    if x.cur_length < 1
        throw(BoundsError())
    end
    return x.values[1]
end

import Base: length
length(x::BoundedPriorityQueue) = x.cur_length

end
