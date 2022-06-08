#Subject
"""
    applylimitrule!(data::PKSubject, rule::LimitRule)

Apply rule to PK subject .

* STEP 1 (NaN step): replace all `NaN` and `missing` values with nan keyword value (if `nan` not NaN);
* STEP 2 (LLOQ step): replace values below `lloq` with `btmax` value if this value befor Tmax or with atmax if this value after Tmax (if `lloq` not NaN);
* STEP 3 (remove NaN): `rm` == true, then remove all `NaN` and `missing` values.
"""
function applylimitrule!(data::PKSubject, rule::LimitRule)
    applylimitrule!(data.time, data.obs, rule)
    data
end
"""
    applylimitrule!(f::Function, data::DataSet{T}, rule::LimitRule) where T <: PKSubject
"""
function applylimitrule!(f::Function, data::DataSet{T}, rule::LimitRule) where T <: PKSubject
    for i in data
        if f(i) applylimitrule!(i, rule) end
    end
    data
end
#DS ind Int
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule, ind::Int) where T <: PKSubject
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule, ind::Int) where T <: PKSubject
    applylimitrule!(data[ind], rule)
    data
end
#DS iter Int
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: PKSubject
    for i in inds
        applylimitrule!(data[i], rule)
    end
    data
end
#DS all
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule) where T <: PKSubject
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule) where T <: PKSubject
    for i = 1:length(data)
        applylimitrule!(data[i], rule)
    end
    data
end
#DS Dict
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule, sort::Dict) where T <: PKSubject
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule, sort::Dict) where T <: PKSubject
    for i = 1:length(data)
        if sort ⊆ data[i].id applylimitrule!(data[i], rule) end
    end
    data
end
"""
    applylimitrule!(time, obs, rule::LimitRule)
"""
function applylimitrule!(time, obs, rule::LimitRule)
    if validobsn(time, obs) == 0 return Float64[], Float64[] end
    cmax, tmax, tmaxn = ctmax(time, obs)
    #NaN Rule
    obsn = length(obs)
    if !isnan(rule.nan)
        for i = 1:obsn
            if isnanormissing(obs[i])
                obs[i] = rule.nan
            end
        end
    end
    #LLOQ rule
    if !isnan(rule.lloq)
        for i = 1:obsn
            if !isnanormissing(obs[i]) && obs[i] <= rule.lloq
                if i <= tmaxn
                    obs[i] = rule.btmax
                else
                    obs[i] = rule.atmax
                end
            end
        end
    end
    #NaN Remove rule
    if rule.rm
        inds = findall(isnanormissing, obs)
        deleteat!(time, inds)
        deleteat!(obs, inds)
    end
    time, obs
end
