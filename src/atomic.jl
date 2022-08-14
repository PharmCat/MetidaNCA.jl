

"""
    cmax(time::AbstractVector, obs::AbstractVector)

Return Cmax
"""
function cmax(time::AbstractVector, obs::AbstractVector)
    length(time) == length(obs) || error("length(time) != length(obs)")
    cmax, tmax, tmaxn = ctmax(time, obs)
    cmax
end

"""
    tmax(time::AbstractVector, obs::AbstractVector)

Return Tmax
"""
function tmax(time::AbstractVector, obs::AbstractVector)
    length(time) == length(obs) || error("length(time) != length(obs)")
    cmax, tmax, tmaxn = ctmax(time, obs)
    tmax
end

"""
    auc(time::AbstractVector, obs::AbstractVector;  calcm = :lint)

Return AUC. All concentration points included in calculation.

* `calcm` - AUC/AUMC calculation method:
    - `:lint` - linear trapezoidal;
    - `:logt` - log-trapezoidal after Tmax;
    - `:luld` - linar up log down;
    - `:luldt` - linear up log down after Tmax;

!!! note
    This function doesn't contain `NaN`, `missing` or dosing time checks.

"""
function auc(time::AbstractVector, obs::AbstractVector;  calcm = :lint)
    length(time) == length(obs) || error("length(time) != length(obs)")
    length(time) >= 2 || error("length(time) >= 2")
    auc_ = 0.0
    if calcm == :lint
        @inbounds for i = 1:(length(time) - 1)
            auc_  += linauc(time[i], time[i + 1], obs[i], obs[i+1])
        end
    else
        cmax, tmax, tmaxn = ctmax(time, obs)
        @inbounds for i = 1:(length(time) - 1)
            auc_  += aucpart(time[i], time[i + 1], obs[i], obs[i + 1], calcm, i >= tmaxn)
        end
    end

    auc_
end
