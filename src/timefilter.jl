
"""
    timefilter(data_::DataSet{<: PKSubject}, time::Tuple{<:Number, <:Number})

Make copy of the data_ and remove all observations < time[1] or > time[2]. Then resize keldata to 0.
"""
function timefilter(data_::DataSet{<: PKSubject}, time::Tuple{<:Number, <:Number})
    data = deepcopy(data_)
    for subj in data
        inds = Int[]
        for n = 1:length(subj)
            if subj.time[n] < time[1] || subj.time[n] > time[2] push!(inds, n) end
        end
        deleteat!(subj.time, inds)
        deleteat!(subj.obs, inds)
        resize!(subj.keldata, 0)
    end
    data
end
