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
    #=
    cmax, tmax, tmaxn = ctmax(data)
    #NaN Rule
    obsn = length(data)
    if !isnan(rule.nan)
        for i = 1:length(data)
            if isnanormissing(data.obs[i])
                data.obs[i] = rule.nan
            end
        end
    end
    #LLOQ rule
    if !isnan(rule.lloq)
        for i = 1:obsn
            if data.obs[i] <= rule.lloq
                if i <= tmaxn
                    data.obs[i] = rule.btmax
                else
                    data.obs[i] = rule.atmax
                end
            end
        end
    end
    #NaN Remove rule
    if rule.rm
        inds = findall(isnanormissing, data.obs)
        deleteat!(data.time, inds)
        deleteat!(data.obs, inds)
    end
    =#
    data
end
#DS ind Int
#DS iter Int
#DS all
#DS Dict



function applylimitrule!(time, obs, rule::LimitRule)
    cmax, tmax, tmaxn = ctmax(time, obs)
    #NaN Rule
    obsn = length(obs)
    if !isnan(rule.nan)
        for i = 1:length(data)
            if isnanormissing(obs[i])
                obs[i] = rule.nan
            end
        end
    end
    #LLOQ rule
    if !isnan(rule.lloq)
        for i = 1:obsn
            if obs[i] <= rule.lloq
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
