using Enzyme: Reverse, Const, Duplicated, autodiff

include("calcium_fluorescence_model.jl")
include("calcium_fluorescence_reparam.jl")

@inline function addGradLogG_fluorescence!(g, k, x_cur, scr)
    @inbounds y = scr.y[k]
    autodiff(Reverse, lG_fluorescence_worker, Const(k), Const(x_cur), Const(y), Const(scr.c), Duplicated(scr.par, g))
    nothing
end

@inline function addGradLogM_calcium!(g, k, x_prev, x_cur, scr)    
    autodiff(Reverse, lM_calcium_worker, Const(k), Const(x_prev), Const(x_cur), Const(scr.c), Duplicated(scr.par, g))
    nothing
end

