export Transformation

abstract Transformation{T,D} <: AbstractSequence{T,D}

for trans in (:Cumprod,:Cumsum,:Aitken,:Euler)
    @eval begin
        export $trans
        type $trans{T,D} <: Transformation{T,D}
            sequence::Vector{T}
            generator::D
        end
        $trans(sequence) = $trans(eltype(sequence)[],sequence)
    end
end

for trans in (:Difference,:Shift)
    @eval begin
        export $trans, order
        type $trans{T,D,O} <: Transformation{T,D}
            sequence::Vector{T}
            generator::D
            function $trans(sequence::Vector{T},generator::D)
                new(sequence,generator)
            end
        end
        $trans(sequence,generator,O::Int) = (@assert O ≥ 0;$trans{eltype(sequence),typeof(generator),O}(sequence,generator))
        $trans(sequence,O::Int) = $trans(eltype(sequence)[],sequence,O)
        $trans(sequence) = $trans(sequence,1)
        order{T,D,O}(::$trans{T,D,O}) = O
        order{T,D,O}(::Type{$trans{T,D,O}}) = O
        *(T1::Type{$trans},T2::$trans) = $trans(sequence(T2),generator(T2),order(T2)+1)
    end
end

export Δ, E
typealias Δ{T,D} Difference{T,D}
typealias E{T,D} Shift{T,D}


for (op,trans) in ((:(Base.cumprod),:Cumprod),(:(Base.cumsum),:Cumsum),(:(Base.diff),:Δ))
    @eval begin
        $op(S::AbstractSequence) = $trans($op(sequence(S)),S)
    end
end

for trans in (:ϵ,:ρ)
    @eval begin
        export $trans
        type $trans{T,D} <: Transformation{T,D}
            sequence::Vector{T}
            generator::D
            auxiliary::Vector{T}
        end
        $trans(sequence) = $trans(eltype(sequence)[],sequence,eltype(sequence)[])
        auxiliary(T::$trans) = T.auxiliary
    end
end

typealias Shanks{T,D} ϵ{T,D}


for trans in (:Drummond,:Levin,:Weniger)
    @eval begin
        export $trans
        type $trans{T,D,Symbol} <: Transformation{T,D}
            sequence::Vector{T}
            generator::D
            numerator::Vector{T}
            denominator::Vector{T}
            function $trans(sequence::Vector{T},generator::D)
                new(sequence,generator,zero(sequence),zero(sequence))
            end
        end
        $trans(sequence,ω::Symbol) = $trans{eltype(sequence),typeof(sequence),ω}(eltype(sequence)[],sequence)
        numerator(T::$trans) = T.numerator
        denominator(T::$trans) = T.denominator
    end
end

for trans in (:Drummond,:Levin)
    @eval begin
        remainder{T,D}(S::$trans{T,D,:u},k::Int) = (s = generator(S);k>1?(k+1)*(s(k)-s(k-1)):2s(1))
        remainder{T,D}(S::$trans{T,D,:t},k::Int) = (s = generator(S);k>1?s(k)-s(k-1):s(1))
        remainder{T,D}(S::$trans{T,D,:d},k::Int) = (s = generator(S);s(k+1)-s(k))
        remainder{T,D}(S::$trans{T,D,:v},k::Int) = (s = generator(S);k>1?(s(k)-s(k-1))*(s(k+1)-s(k))/(2s(k)-s(k-1)-s(k+1)):s(1)*(s(2)-s(1))/(2s(1)-s(2)))
    end
end

remainder{T,D}(S::Weniger{T,D,:y},k::Int) = (s = generator(S);k>1?(k+1)*(s(k)-s(k-1)):2s(1))
remainder{T,D}(S::Weniger{T,D,:τ},k::Int) = (s = generator(S);k>1?s(k)-s(k-1):s(1))
remainder{T,D}(S::Weniger{T,D,:δ},k::Int) = (s = generator(S);s(k+1)-s(k))
remainder{T,D}(S::Weniger{T,D,:φ},k::Int) = (s = generator(S);k>1?(s(k)-s(k-1))*(s(k+1)-s(k))/(2s(k)-s(k-1)-s(k+1)):s(1)*(s(2)-s(1))/(2s(1)-s(2)))

Drummond(sequence) = Drummond(sequence,:u)
Levin(sequence) = Levin(sequence,:u)
Weniger(sequence) = Weniger(sequence,:y)
