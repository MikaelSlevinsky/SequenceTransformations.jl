module SequenceTransformations
    using Base
    import Base: +,-,.+,.-,*,.*,/,./

include("Sequence.jl")
include("Transformations.jl")
include("generate.jl")

end # module
