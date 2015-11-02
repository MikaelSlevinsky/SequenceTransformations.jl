# Generating new terms in a sequence

function generate!(S::Sequence,n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n s[k] = a(k) end
end

function generate!(S::Cumsum,n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n s[k] = a(k)+s[k-1] end
end

function generate!(S::Cumprod,n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n s[k] = a(k)*s[k-1] end
end

function generate!{T,D}(S::Difference{T,D,1},n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n s[k] = a(k+1)-a(k) end
end

function generate!{T,D}(S::Difference{T,D,2},n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n s[k] = a(k+2)-2a(k+1)+a(k) end
end

function generate!{T,D,O}(S::Shift{T,D,O},n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n s[k] = a(k+O) end
end


# Aitken's Δ² process

function generate!(S::Aitken,n::Int)
    s,a = sequence(S),generator(S)
    i = length(s)
    i ≤ n && resize!(s,n)
    for k = 1+i:n
        s[k] = a(k+2)-(a(k+2)-a(k+1))^2/(a(k+2)-2a(k+1)+a(k))
    end
end


# Wynn's ϵ- and ρ-algorithms

for trans in (:ϵ,:ρ)
    rec = parse(string(trans)*"recurrence!")
    @eval begin
        function generate!(S::$trans,n::Int)
            data,s = sequence(S),generator(S)
            aux = auxiliary(S)
            i = length(data)
            i ≤ n && resize!(data,n);resize!(aux,n)
            for k = 1+i:n
                aux[k] = s(k)
                $rec(k,aux,data)
            end
        end
    end
end

function ϵrecurrence!{T}(k::Int,aux::Vector{T},data::Vector{T})
    aux2 = zero(T)
    for j = k:-1:2
        aux1 = aux2
        aux2 = aux[j-1]
        diff = aux[j]-aux2
        aux[j-1] = aux1 + one(T)/diff
    end
    data[k] = isodd(k) ? aux[1] : aux[2]
end

function ρrecurrence!{T}(k::Int,aux::Vector{T},data::Vector{T})
    aux2 = zero(T)
    for j = k:-1:2
        aux1 = aux2
        aux2 = aux[j-1]
        diff = aux[j]-aux2
        aux[j-1] = aux1 + j/diff
    end
    data[k] = isodd(k) ? aux[1] : aux[2]
end


# Drummond, Levin, and Weniger transformations


for trans in (:Drummond,:Levin,:Weniger)
    rec = parse(string(trans)*"recurrence!")
    @eval begin
        function generate!{T}(S::$trans{T},n::Int)
            ℓ,s = sequence(S),generator(S)
            N,D = numerator(S),denominator(S)
            i = length(ℓ)
            i ≤ n && resize!(ℓ,n);resize!(N,n);resize!(D,n)
            for k = 1+i:n
                ω = remainder(S,k)
                N[k],D[k] = s(k)/ω,one(T)/ω
                $rec(k,ω,N,D,ℓ)
            end
        end
    end
end

function Drummondrecurrence!{T}(k::Int,ω::T,N::Vector{T},D::Vector{T},ℓ::Vector{T})
    for j = 1:k-1
        N[k-j] = N[k-j+1]-N[k-j]
        D[k-j] = D[k-j+1]-D[k-j]
    end
    ℓ[k] = N[1]/D[1]
end

function Levinrecurrence!{T}(k::Int,ω::T,N::Vector{T},D::Vector{T},ℓ::Vector{T})
    kp1 = k+one(T)
    kdkp1 = k/kp1
    kdkp1jm2 = 1/kdkp1
    for j = 1:k-1
        cst = -(kp1-j)*kdkp1jm2/kp1
        N[k-j] = muladd(cst,N[k-j],N[k-j+1])
        D[k-j] = muladd(cst,D[k-j],D[k-j+1])
        kdkp1jm2 *= kdkp1
    end
    ℓ[k] = N[1]/D[1]
end

function Wenigerrecurrence!{T}(k::Int,ω::T,N::Vector{T},D::Vector{T},ℓ::Vector{T})
    km1 = k-one(T)
    for j = 1:k-1
        kjm1 = km1+j
        kjm2 = kjm1-1
        cst = -k*(km1/kjm1)/kjm2
        N[k-j] = muladd(cst,N[k-j],N[k-j+1])
        D[k-j] = muladd(cst,D[k-j],D[k-j+1])
    end
    ℓ[k] = N[1]/D[1]
end
