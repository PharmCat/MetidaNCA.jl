function dropfunction!(f::Function, subj::PKSubject)
    inds = findall(f, subj.obs)
    deleteat!(subj.time, inds)
    deleteat!(subj.obs, inds)
    return subj
end

function dropfunction(f::Function, subj::PKSubject)
    return dropfunction!(f, deepcopy(subj))
end

dropnanormissinf!(subj::PKSubject)    =  dropfunction!(isnanormissing, subj)
dropnanormissinf(subj::PKSubject)     =  dropfunction(isnanormissing, subj)

dropnan!(subj::PKSubject)    =  dropfunction!(isnan, subj)
dropnan(subj::PKSubject)     =  dropfunction(isnan, subj)

dropmissing!(subj::PKSubject) =  dropfunction!(ismissing, subj)
dropmissing(subj::PKSubject)  =  dropfunction(ismissing, subj)



