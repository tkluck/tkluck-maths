module Heaps

#heapparent(ix) = div(ix, 2)
#heapleft(ix) = 2*ix
#heapright(ix) = heapleft(ix) + 1

const ChunkSize = 8 # 2^16
const ChunkFanout = div(ChunkSize, 2)

function heapparent(ix::Int)
    chunk = div(ix, ChunkSize)
    if chunk == 0
        return div(ix, 2)
    end
    chunk_base0 = chunk*ChunkSize - 1
    ix_in_chunk = ix - chunk_base0

    if ix_in_chunk > 4
        return chunk_base0 + div(ix_in_chunk-1, 2)+1
    elseif ix_in_chunk > 2
        return ix - 2
    else
        ChunkFanout = div(ChunkSize+1, 2)
        parent_chunk = div(chunk-1, ChunkFanout)
        chunk_index_as_child = chunk - ChunkFanout * parent_chunk
        n = (parent_chunk+1) * ChunkSize - ChunkFanout + chunk_index_as_child - 1
        return n
    end
end

function heapleft(ix::Int)
    chunk = div(ix, ChunkSize)
    if ix < ChunkFanout
        return 2*ix
    elseif chunk == 0
        child_chunk_ix = ix - ChunkFanout + 1
        return child_chunk_ix*ChunkSize
    end

    chunk_base0 = chunk*ChunkSize - 1
    ix_in_chunk = ix - chunk_base0

    if ix_in_chunk <= 2
        return ix + 2
    elseif ix_in_chunk <= ChunkFanout
        return chunk_base0 + 2*(ix_in_chunk-1) + 1
    else
        ix_as_child = ix_in_chunk - ChunkFanout
        child_chunk_ix = ChunkFanout*chunk + ix_as_child
        return child_chunk_ix*ChunkSize
    end

end
function heapright(ix::Int)
    if ix >= ChunkSize && (ix % ChunkSize) <= 1
        return heapleft(ix)
    else
        return heapleft(ix) + 1
    end
end

function percolate_up!(xs::AbstractArray, i::Integer, o::Base.Ordering)
    x = xs[i]
    @inbounds while (j = heapparent(i)) >= 1
        if Base.Order.lt(o, x, xs[j])
            xs[i] = xs[j]
            i = j
        else
            break
        end
    end
    xs[i] = x
end

function percolate_down!(xs::AbstractArray, i::Integer, o::Base.Ordering, len::Integer)
    x = xs[i]
    @inbounds while (l = heapleft(i)) <= len
        r = heapright(i)
        j = r > len || Base.Order.lt(o, xs[l], xs[r]) ? l : r
        if Base.Order.lt(o, xs[j], x)
            xs[i] = xs[j]
            i = j
        else
            break
        end
    end
    xs[i] = x
end

function sort!(xs::AbstractArray)
    for i = heapparent(length(xs)):-1:1
        percolate_down!(xs, i, Base.Order.Reverse, length(xs))
    end
    for i = length(xs):-1:2
        xs[1], xs[i] = xs[i], xs[1]
        percolate_down!(xs, 1, Base.Order.Reverse, i-1)
    end
end

end
