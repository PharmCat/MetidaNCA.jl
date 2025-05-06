
"""
    timefilter(subj::PKSubject, time::AbstractRange)

Exclude all observation than not in time range.
"""
function timefilter(subj::PKSubject, time::AbstractRange)
    subj_ = deepcopy(subj)
    inds = Int[]
    for n = 1:length(subj_)
        if !(gettime(subj_)[n] in time) push!(inds, n) end
    end
    deleteat!(gettime(subj_), inds)
    deleteat!(getobs(subj_), inds)
    resize!(subj_.keldata, 0)
    if !(subj_.kelrange.kelstart in time) || !(subj_.kelrange.kelend in time) || any(x-> !(x in time), subj_.kelrange.kelexcl)
        subj_.kelrange = ElimRange()
        subj_.kelauto = true
    end
    subj_
end

"""
    timefilter(subj::PKSubject, time::Tuple{<:Number, <:Number})

Make deepcopy of subj and remove all observations < time[1] or > time[2]. Then resize keldata to 0.

If any of points in elimination rage not in min/max time, then elimination settings reset.
"""
function timefilter(subj::PKSubject, time::Tuple{<:Number, <:Number})
    timefilter(subj, LinRange(time[1], time[2], 2))
end
"""
    timefilter(data::DataSet{<: PKSubject}, time)

Make new DataSet with new filtered subjects.
"""
function timefilter(data::DataSet{<: PKSubject}, time)
    subj  = getdata(data)
    data_ = similar(subj)
    for i in 1:length(subj)
        data_[i] = timefilter(subj[i], time)
    end
    DataSet(data_)
end
