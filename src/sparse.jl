"""
    auc_sparse(time, obs)

AUC for sparse data.

```math
w_1 = (t_2 - t_1) / 2
```

```math
w_j = (t_{j+1} - t_{j-1}) / 2  (2 \\leq j \\leq J - 1)
```

```math
w_J = (t_J - t_{J-1}) / 2
```

```math
AUC = \\sum_{j=1}^J \\mu_j w_j
```

where `math \\mu_j` is the mean drug concentration at time t.
"""
function auc_sparse(time::AbstractVector, obs::AbstractVector)
    if length(time) < 2 error("length(time) < 2") end
    if length(time) != length(obs) error("length(time) != length(obs)") end
    for i = 1:length(time) - 1
        if time[i+1] <= time[i] error("Unsorted observations!") end
    end
    wts = Vector{Float64}(undef, length(time))
    wts[1]   = (time[2] - time[1]) / 2
    wts[end] = (time[end] - time[end-1]) / 2
    if length(time) > 2
        for i = 2:length(time) - 1
            wts[i] = (time[i+1] - time[i-1]) / 2
        end
    end
    return obs' * wts
end

