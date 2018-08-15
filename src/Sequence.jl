export AbstractSequence, Sequence, sequence, generator

abstract type AbstractSequence{T,D} end

struct Generator{D}
    generator::D
end

(G::Generator)(i) = G.generator(i)
Base.getindex(G::Generator, i) = G(i)

mutable struct Sequence{T,D} <: AbstractSequence{T,D}
    sequence::Vector{T}
    generator::Generator{D}
end

Sequence(generator::D) where D<:Base.Callable = Sequence([generator(1)],Generator(generator))

sequence(S::AbstractSequence) = S.sequence
generator(S::AbstractSequence) = S.generator

Base.eltype(::AbstractSequence{T}) where T = T
Base.size(::AbstractSequence) = (Inf,)
Base.length(::AbstractSequence) = Inf

Base.getindex(S::AbstractSequence,i::Integer) = (generate!(S,i); sequence(S)[i])
Base.getindex(S::AbstractSequence,ir::AbstractRange) = (generate!(S,maximum(ir)); sequence(S)[ir])
