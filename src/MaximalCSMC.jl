module MaximalCSMC

using RNGPool
import Statistics.mean
using Random

include("structures.jl")
include("subroutines.jl")
include("algorithms.jl")
include("unbiased.jl")
include("couplingTime.jl")
include("optimizers.jl")
include("gradientDescent.jl")

# The specific forward coupling strategies are implemented in these:
include("forwardCouplings/filterStateMaximalCoupling.jl")
include("forwardCouplings/particleMaximalCoupling.jl")
include("forwardCouplings/indexMaximalCoupling.jl")
include("forwardCouplings/jointIndexMaximalCoupling.jl")

export CCSMCIO, ccXpf!, ccMCBSpf!, initializeCCSMC, SGML, FilterStateMaximalCoupling, 
  ParticleMaximalCoupling, IndexMaximalCoupling, HybridMaximalCoupling, HybridLagMaximalCoupling,
  JointIndexCoupling

end # module
