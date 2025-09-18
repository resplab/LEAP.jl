
PROVINCE_MAP = Dict(
    "Canada" => "CA",
    "British Columbia" => "BC",
    "Alberta" => "AB",
    "Saskatchewan" => "SK",
    "Manitoba" => "MB",
    "Ontario" => "ON",
    "Quebec" => "QC",
    "Newfoundland and Labrador" => "NL",
    "Nova Scotia" => "NS",
    "New Brunswick" => "NB",
    "Prince Edward Island" => "PE",
    "Yukon" => "YT",
    "Northwest Territories" => "NT",
    "Nunavut" => "NU"
)

function get_province_id(province)
    return PROVINCE_MAP[province]
end


function get_sex_id(sex)
    return sex[1:1]
end


function format_age_group(age_group::AbstractString, upper_age_group="100 years and over")
    if age_group == upper_age_group
        age = replace(age_group, " years and over" => "")
        age = tryparse(Int, age)
    elseif age_group == "Under 1 year"
        age = 0
    else
        age = replace(age_group, " years" => "")
        age = replace(age, " year" => "")
        age = tryparse(Int, age)
    end
    return age
end
