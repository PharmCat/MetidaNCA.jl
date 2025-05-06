function dropfunction!(f::Function, subj::PKSubject)
    inds = findall(f, getobs(subj))
    deleteat!(gettime(subj), inds)
    deleteat!(getobs(subj), inds)
    return subj
end

function dropfunction(f::Function, subj::PKSubject)
    return dropfunction!(f, deepcopy(subj))
end

dropnanormissing!(subj::PKSubject)    =  dropfunction!(isnanormissing, subj)
dropnanormissing(subj::PKSubject)     =  dropfunction(isnanormissing, subj)

dropnan!(subj::PKSubject)    =  dropfunction!(isnan, subj)
dropnan(subj::PKSubject)     =  dropfunction(isnan, subj)

dropmissing!(subj::PKSubject) =  dropfunction!(ismissing, subj)
dropmissing(subj::PKSubject)  =  dropfunction(ismissing, subj)



