export AbstractSequence, Sequence, sequence, generator

abstract AbstractSequence{T,D}

type Sequence{T,D} <: AbstractSequence{T,D}
    sequence::Vector{T}
    generator::D
end

Sequence{D<:Base.Callable}(generator::D) = Sequence([generator(1)],generator)

sequence(S::AbstractSequence) = S.sequence
generator(S::AbstractSequence) = S.generator

Base.eltype{T}(::AbstractSequence{T}) = T
Base.size(::AbstractSequence) = (Inf,)
Base.length(::AbstractSequence) = Inf
Base.start{T,D}(S::AbstractSequence{T,D}) = S(1)
Base.next{T,D}(S::AbstractSequence{T,D},i::Int) = S(i),i+1
Base.next{T,D}(S::AbstractSequence{T,D},i::AbstractFloat) = next(S,round(Int,i))
Base.done{T,D}(S::AbstractSequence{T,D},i) = i > length(sequence(S))

Base.call{T,D}(S::AbstractSequence{T,D},i::Int) = (generate!(S,i);sequence(S)[i])
Base.call{T,D}(S::AbstractSequence{T,D},ir::Range) = (generate!(S,maximum(ir));sequence(S)[ir])
Base.getindex(S::AbstractSequence,i) = S(i)
