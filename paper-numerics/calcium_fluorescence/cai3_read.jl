using MAT

function cai3_read(data_id, cell_id; prefix="", drop_last=1)
    vars = matread("$(prefix)data.$(data_id).train.preprocessed.mat")
    data = vars["data"][cell_id]
    
    F = data["calcium"][1:end-drop_last]
        
    # Scale to [0,1]:
    F = F .- minimum(F)
    F ./= maximum(F)

    N = Int.(data["spikes"][1:end-drop_last])

    # Time step:
    Δ = 1/data["fps"]

    # Time of data
    T = length(F)
    Ts = (1:T)*Δ

    F, N, Δ, Ts
end


