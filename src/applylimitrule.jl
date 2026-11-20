#Subject
"""
    applylimitrule!(data::Union{PKSubject, PDSubject}, rule::LimitRule)

Apply rule to PK subject .

* STEP 1: if `replace_blq` NaN or values < lloq replaced btmax and atmax values

* STEP 2: if `drop_missing` drop all missing values;
"""
function applylimitrule!(data::Union{PKSubject, PDSubject}, rule::LimitRule)
    if rule.drop_missing && length(getfield(data, :obs)) != 1 error("Can't remove observations if multiple data provided.") end
    applylimitrule!(data.time, getobs(data), rule)
    data
end
"""
    applylimitrule!(f::Function, data::DataSet{T}, rule::LimitRule) where T <: Union{PKSubject, PDSubject}

Apply if `f(subj)` return  `true`.
"""
function applylimitrule!(f::Function, data::DataSet{T}, rule::LimitRule) where T <: Union{PKSubject, PDSubject}
    for i in data
        if f(i) applylimitrule!(i, rule) end
    end
    data
end
#DS ind Int
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule, ind::Int) where T <: Union{PKSubject, PDSubject}

Apply by ind.
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule, ind::Int) where T <: Union{PKSubject, PDSubject}
    applylimitrule!(data[ind], rule)
    data
end
#DS iter Int
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: Union{PKSubject, PDSubject}

Apply by inds.
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule, inds::Union{Vector{Int}, UnitRange{Int}, Tuple{Vararg{Int}}}) where T <: Union{PKSubject, PDSubject}
    for i in inds
        applylimitrule!(data[i], rule)
    end
    data
end
#DS all
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule) where T <: Union{PKSubject, PDSubject}

Apply to all dataset.
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule) where T <: Union{PKSubject, PDSubject}
    for i = 1:length(data)
        applylimitrule!(data[i], rule)
    end
    data
end
#DS Dict
"""
    applylimitrule!(data::DataSet{T}, rule::LimitRule, sort::Dict) where T <: PKSubject
"""
function applylimitrule!(data::DataSet{T}, rule::LimitRule, sort::Dict) where T <: Union{PKSubject, PDSubject}
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
    #Replace BLQ  rule
    if !isnan(rule.replace_blq)
        for i = 1:obsn
            if !ismissing(obs[i]) && (isnan(obs[i]) || obs[i] <= rule.lloq)
                if i <= tmaxn
                    obs[i] = rule.btmax
                else
                    obs[i] = rule.atmax
                end
            end
        end
    end
    #Missing Remove rule
    if rule.drop_missing
        inds = findall(ismissing, obs)
        deleteat!(time, inds)
        deleteat!(obs, inds)
    end
    time, obs
end
