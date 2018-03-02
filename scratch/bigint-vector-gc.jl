function busywork(container,n=1000)
    for i=1:n
        container = container .+ container[end:-1:1]
    end
    container
end

version() = VersionNumber(unsafe_string(unsafe_load(cglobal((:__gmp_version, :libgmp), Ptr{Cchar}))))
bits_per_limb() = Int(unsafe_load(cglobal((:__gmp_bits_per_limb, :libgmp), Cint)))

const VERSION = version()
const BITS_PER_LIMB = bits_per_limb()

# GMP's mp_limb_t is by default a typedef of `unsigned long`, but can also be configured to be either
# `unsigned int` or `unsigned long long int`. The correct unsigned type is here named Limb, and must
# be used whenever mp_limb_t is in the signature of ccall'ed GMP functions.
if BITS_PER_LIMB == 32
    const Limb = UInt32
    const SLimbMax = Union{Int8, Int16, Int32}
    const ULimbMax = Union{UInt8, UInt16, UInt32}
elseif BITS_PER_LIMB == 64
    const Limb = UInt64
    const SLimbMax = Union{Int8, Int16, Int32, Int64}
    const ULimbMax = Union{UInt8, UInt16, UInt32, UInt64}
else
    error("GMP: cannot determine the type mp_limb_t (__gmp_bits_per_limb == $BITS_PER_LIMB)")
end


mutable struct MyBigInt <: Signed
    alloc::Cint
    size::Cint
    d::Ptr{Limb}


    function MyBigInt()
        b = new(zero(Cint), zero(Cint), C_NULL)
        ccall((:__gmpz_init, :libgmp), Void, (Ref{MyBigInt},), b)
        return b
    end

    function MyBigInt(a)
        b = new(zero(Cint), zero(Cint), C_NULL)
        ccall((:__gmpz_init, :libgmp), Void, (Ref{MyBigInt},), b)
        ccall((:__gmpz_set, :libgmp), Void, (Ref{MyBigInt}, Ref{BigInt}), b, BigInt(a))
        return b
    end
end

MyBigInt(a::MyBigInt) = a

import Base: +
function +(a::MyBigInt, b::MyBigInt)
    res = MyBigInt()
    ccall((:__gmpz_add, :libgmp), Void, (Ref{MyBigInt}, Ref{MyBigInt}, Ref{MyBigInt}), res, a, b)
    return res
end

const mpz_t = Ref{BigInt}
const my_mpz_t = Ref{MyBigInt}

mutable struct BigIntVector{V}
    x::V
end

function BigIntVector(iter)
    x = map(MyBigInt, iter)
    v = BigIntVector{typeof(x)}(x)

    finalizer(v) do v
        for b in v.x
            ccall((:__gmpz_clear, :libgmp), Void, (Ref{MyBigInt},), b)
        end
    end
    v
end


import Base: lastindex, getindex, broadcast

lastindex(v::BigIntVector) = lastindex(v.x)
getindex(v::BigIntVector, ix...) = getindex(v.x, ix...)
broadcast(f, v::BigIntVector, args...) = typeof(v)(broadcast(f, v.x, args...)

busywork(collect(1:100))
busywork(BigInt.(1:100))
busywork(BigIntVector(BigInt.(1:100)))

@time busywork(collect(1:10000))
@time busywork(BigInt.(1:10000))
@time busywork(BigIntVector(BigInt.(1:10000)))

