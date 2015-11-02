export Transformation, LinearTransformation, NonlinearTransformation

abstract Transformation{T,D} <: AbstractSequence{T,D}

abstract LinearTransformation{T,D} <: Transformation{T,D}

abstract NonlinearTransformation{T,D} <: Transformation{T,D}

for trans in (:Cumprod,:Cumsum)
    @eval begin
        export $trans
        type $trans{T,D} <: LinearTransformation{T,D}
            sequence::Vector{T}
            generator::D
        end
        $trans(generator) = $trans(eltype(generator)[],generator)
    end
end

for trans in (:Difference,:Shift)
    @eval begin
        export $trans, order
        type $trans{T,D,O} <: LinearTransformation{T,D}
            sequence::Vector{T}
            generator::D
            function $trans(sequence::Vector{T},generator::D)
                new(sequence,generator)
            end
        end
        $trans(sequence,generator,O::Int) = (@assert O ≥ 0;$trans{eltype(sequence),typeof(generator),O}(sequence,generator))
        $trans(generator,O::Int) = $trans(eltype(generator)[],generator,O)
        $trans(generator) = $trans(generator,1)
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

for trans in (:Aitken,:Euler)
    @eval begin
        export $trans
        type $trans{T,D} <: NonlinearTransformation{T,D}
            sequence::Vector{T}
            generator::D
        end
        $trans(generator) = (D = eltype(generator); T = typeof(one(D)/one(D));$trans(T[],generator))
    end
end

for trans in (:ϵ,:ρ)
    @eval begin
        export $trans
        type $trans{T,D} <: NonlinearTransformation{T,D}
            sequence::Vector{T}
            generator::D
            auxiliary::Vector{T}
        end
        $trans(generator) = (D = eltype(generator); T = typeof(one(D)/one(D));$trans(T[],generator,T[]))
        auxiliary(T::$trans) = T.auxiliary
    end
end

typealias Shanks{T,D} ϵ{T,D}


for trans in (:Drummond,:Levin,:Weniger)
    @eval begin
        export $trans
        type $trans{T,D,Symbol} <: NonlinearTransformation{T,D}
            sequence::Vector{T}
            generator::D
            numerator::Vector{T}
            denominator::Vector{T}
            function $trans(sequence::Vector{T},generator::D)
                new(sequence,generator,zero(sequence),zero(sequence))
            end
        end
        $trans(generator,ω::Symbol) = (D = eltype(generator);T = typeof(one(D)/one(D)); $trans{T,typeof(generator),ω}(T[],generator))
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

Drummond(generator) = Drummond(generator,:u)
Levin(generator) = Levin(generator,:u)
Weniger(generator) = Weniger(generator,:y)
