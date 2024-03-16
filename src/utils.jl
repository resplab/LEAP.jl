function dict_initializer(parameter_names::Union{Nothing,Vector{Symbol}})
    isnothing(parameter_names) ? nothing : Dict(parameter_names .=> missing)
end

function vec_to_dict(v::AbstractArray, ll::AbstractVector)::AbstractDict
    d = Dict()
        for i in eachindex(ll)
            d[ll[i]] = v[i]
        end
    return d
end


function string_to_symbols_dict(dict::AbstractDict)::AbstractDict
    new_dict = Dict()
    for (key, value) in dict
        new_dict[Symbol(key)] = value
    end
    return new_dict
end
