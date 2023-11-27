"""
    Location

A struct containing the location information of a person (agent).

# Fields
- `census_division::String`: the name of the census division to which the person belongs.
- `longitude::Union{Nothing, Float64}`: the longitude of the person's home.
- `latitude::Union{Nothing, Float64}`: the latitude of the person's home.
"""
struct Location
    census_division::Union{Nothing, String}
    longitude::Union{Nothing, Float64}
    latitude::Union{Nothing, Float64}
end


function process_initial(agent::Agent)
    return
end
