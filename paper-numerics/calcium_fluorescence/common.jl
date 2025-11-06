using ProgressMeter

function aggregate_window(fun, x, n)
    n_x = length(x)
    y = similar(x)
    for i = 1:n_x
        i_start = max(1, i-n)
        i_end = min(n_x, n+i)
        y[i] = fun(x[i_start:i_end])
    end
    y
end

moving_average(x, n) = aggregate_window(mean, x, n)
moving_sum(x, n) = aggregate_window(sum, x, n)

function run_cbpf(model, smcio, lM; n_CBPF = 500, burnin_CBPF = 0.1n_CBPF, showprogress::Bool=true)
    X_ref  = [CalciumParticle() for t in 1:T]
    X_ref_  = [CalciumParticle() for t in 1:T]
    X_smc = []; accepted = zeros(Int64, T)
    smc!(model, smcio)
    SequentialMonteCarlo.pickParticle!(X_ref_, smcio)
    @showprogress showprogress for i = 1:n_CBPF
        csmc!(model, smcio, X_ref_, X_ref)
        SequentialMonteCarlo.pickParticleBS!(X_ref, smcio, lM)
        i > burnin_CBPF && push!(X_smc, deepcopy(X_ref))
        accepted .+= X_ref .!= X_ref_
        X_ref_, X_ref = X_ref, X_ref_
    end
    X_smc, accepted
end

