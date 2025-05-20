using RNGPool
using SMCExamples.LinearGaussian
import SMCExamples.LinearGaussian: LGTheta, simulateLGModel, kalman,
  makeLGModel
import SMCExamples.Particles.Float64Particle

function setupLGModel(n::Int64, ar::Float64, obsVariance::Float64 = 1.0, mutationVariance::Float64 = 1.0)
  statVar = mutationVariance/(1-ar^2)
  theta = LGTheta(ar, mutationVariance, 1.0, obsVariance, 0.0, statVar)
  ys = simulateLGModel(theta, n)
  model = makeLGModel(theta, ys)
  lM = LinearGaussian.makelM(theta)
  ko = kalman(theta, ys)
  return model, lM, ko
end

function setupLG0Model(n::Int64, ar::Float64, obsVariance::Float64=1.0, mutationVariance::Float64 = 1.0)
  statVar = mutationVariance/(1-ar^2)
  theta = LGTheta(ar, mutationVariance, 1.0, obsVariance, 0.0, statVar)
  ys = Vector{Float64}(undef, n)
  for i in 1:n
    ys[i] = 0.0
  end
  model = makeLGModel(theta, ys)
  lM = LinearGaussian.makelM(theta)
  ko = kalman(theta, ys)
  return model, lM, ko
end
