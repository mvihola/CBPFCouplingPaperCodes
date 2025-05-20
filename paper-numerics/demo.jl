using MaximalCSMC

function countEqual(ps1::Vector{Particle},ps2::Vector{Particle}) where Particle
  N = length(ps1)
  count = 0
  for i in 1:N
    count += ps1[i] == ps2[i]
  end
  return count
end

function goodSetProportion(ccsmcio)
  ps1 = ccsmcio.smcio1.allZetas
  ps2 = ccsmcio.smcio2.allZetas
  n = ccsmcio.n
  N = ccsmcio.N
  good = Vector{Float64}(undef, n)
  for i in 1:n
    good[i] = countEqual(ps1[i],ps2[i])/N
  end
  return good
end

function makeMatrix(vs::Vector{Vector{T}}) where T
  m = length(vs)
  n = length(vs[1])
  M = Matrix{T}(undef, m, n)
  for i in 1:m
    M[i,:] = vs[i]
  end
  return M
end

function makeVector(v1::Vector{Particle}, v2::Vector{Particle}) where Particle
  n = length(v1)
  v = Vector{Bool}(undef, n)
  for j in 1:n
    v[j] = v1[j] == v2[j]
  end
  return v
end

@inline function makeCharacter(v1::Vector{Particle}, v2::Vector{Particle},
  start::Int64, finish::Int64) where Particle
  for j in start:finish
    if v1[j] != v2[j]
      return "x"
    end
  end
  return "-"
end

function makeString(v1::Vector{Particle}, v2::Vector{Particle},
  quantum::Int64) where Particle
  s::String = ""
  n::Int64 = length(v1)
  for i = 1:n
    if mod(i, quantum) == 1 || quantum == 1
      start = i
      finish = min(i+quantum-1, n)
      s *= makeCharacter(v1, v2, start, finish)
    end
  end
  s *= "\n"
  return s
end

# run the ccMCBSpf and visualize until the meeting time
function visualizeCCSMCgen(model, lM::F, N::Int64;
  forwardCoupling=ParticleMaximalCoupling(), 
  independentInitialization=false,
  maxit::Int64 = typemax(Int64)) where F<:Function

  println("\nForward with $(typeof(forwardCoupling)), backward index coupling, N = ", N, ":\n")

  ccsmcio = CCSMCIO{model.particle, model.pScratch}(N, model.maxn)

  initializeCCSMC(model, lM, ccsmcio, independentInitialization)

  quantum = max(1, ceil(Int64, ccsmcio.n / displaysize(stdout)[2]))

  ref1 = ccsmcio.ref1
  ref2 = ccsmcio.ref2

  printstyled(makeString(ref1, ref2, quantum), color=:green)

  vs = Vector{Vector{Bool}}(undef, 0)
  push!(vs, makeVector(ref1, ref2))
  props = Vector{Vector{Float64}}(undef, 0)
  push!(props, goodSetProportion(ccsmcio))

  for i in 1:maxit
    ccMCBSpf!(model, ccsmcio, lM; forwardCoupling=forwardCoupling)
    printstyled(makeString(ref1, ref2, quantum), color=:green)
    push!(vs, makeVector(ref1, ref2))
    push!(props, goodSetProportion(ccsmcio))
    if MaximalCSMC.checkEqual(ref1, ref2)
      println("coupled at iteration $i")
      return makeMatrix(vs), makeMatrix(props)
    end
  end
  println("never coupled in $maxit iterations")
  return makeMatrix(vs)
end
