# Implementation of the Adam minimizer (Algorithm 1 of https://arxiv.org/abs/1412.6980)
struct ADAM{T,FT}
    θ::T
    m::T
    v::T
    ε::FT
    α::FT
    β_1::FT
    β_2::FT
end

"""
   ADAM(θ_0; α=0.001, ε=1e-8, β_1=0.9, β_2=0.999)

Initialize Adam minimizer with initial value θ_0 (mutable vector) and optionally parameter values.

# Example

```
θ_0 = randn(100)*100   # Random initial value
o = ADAM(θ_0)          # Initialize the minimizer
g = zeros(100)         # Storage for stochastic gradient
for t = 1:1_000_000
   randn!(g)           # Randomness
   g .+= o.θ           # ...and true gradient value
   step!(o, g, t)      # Step of the Adam minimizer
end
maximum(abs.(o.θ))     # Check how far we are
```
"""
function ADAM(θ_0; α=0.001, ε=1e-8, β_1=0.9, β_2=0.999)
    FT = eltype(θ_0)
    m = similar(θ_0); m .= 0.0
    v = similar(θ_0); v .= 0.0
    ADAM(copy(θ_0), m, v, FT(ε), FT(α), FT(β_1), FT(β_2))
end

function step!(o::ADAM{T,FT}, g::AbstractVector, t::Real) where {T,FT}
    o.m .*= o.β_1 
    o.m .+= (one(FT) - o.β_1) .* g
    o.v .*= o.β_2 
    o.v .+= (one(FT) - o.β_2) .* g.^2
    o.θ .-= o.α .* o.m ./ (one(FT) - o.β_1^t) ./ (sqrt.(o.v ./ (one(FT) - o.β_2^t)) .+ o.ε)
    nothing
end

# Stochastic gradient descent
struct SGD{T,FT}
    θ::T
    c::FT
    γ::FT
end
function SGD(θ_0; c=1.0, γ=0.66)
    FT = eltype(θ_0)
    SGD(θ_0, FT(c), FT(γ))
end

function step!(o::SGD{T,FT}, g::AbstractVector, t::Real) where {T,FT}
    η = o.c * t^(-o.γ)
    o.θ .-= η * g
    nothing
end
