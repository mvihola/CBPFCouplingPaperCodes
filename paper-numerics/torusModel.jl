using SequentialMonteCarlo
using RNGPool
import SMCExamples.Particles.Float64Particle

function setupTorusModel(n::Int64, w::Float64, a::Float64, b::Float64)

  @inline function M!(newParticle::Float64Particle, rng::RNG, p::Int64, particle::Float64Particle, ::Nothing)
    if p == 1
      # Start with uniform
      newParticle.x = rand(rng)
    else
      if rand(rng) > a # With probability 1-a, U(x-w/2, x+w/2) (mod 1)
        newParticle.x = mod(particle.x + w*(rand(rng)-0.5), 1.0)
      else # otherwise, U(0,1)
        newParticle.x = rand(rng)
      end
    end
  end

  @inline function lM(::Int64, particle::Float64Particle,
    newParticle::Float64Particle, ::Nothing)
    diff::Float64 = abs(particle.x - newParticle.x)
    if diff > 0.5
      diff = 1.0 - diff
    end
    if diff <= w/2
      return log(a + (1.0 - a)/w)
    else
      return log(a)
    end
  end
  @inline function lG(::Int64, particle::Float64Particle, ::Nothing)
    if particle.x <= 0.25 || 0.5 < particle.x <= 0.75
      return log(b)
    else # 0.25 < particle.x <= 0.5 || 0.75 < particle.x
      return log(1.0 - b)
    end
  end
  
  return SMCModel(M!, lG, n, Float64Particle, Nothing), lM
end
