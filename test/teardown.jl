using Coverage, LEAP
using ArgParse


function process_coverage_files(module_name="all")
    LEAP_SRC_PATH = dirname(pathof(LEAP))
    TEST_PATH = @__DIR__
    # process '*.cov' files
    if module_name == "all"
        folder = LEAP_SRC_PATH
        coverage = process_folder(folder)
        files = readdir(folder)
        files = map!(file -> joinpath(folder, file), files, files)
        append!(files, readdir(joinpath(folder, "modules")))
        files = map!(file -> joinpath(folder, "modules", file), files, files)
        filter!(file -> occursin(".cov", file), files)
        sources = Vector{String}()
        for file in files
            push!(sources, read(file, String))
        end
    else
        if module_name * ".jl" in readdir(joinpath(LEAP_SRC_PATH, "modules"))
            folder = joinpath(LEAP_SRC_PATH, "modules")
        else
            folder = LEAP_SRC_PATH
        end
        coverage = process_file(joinpath(folder, module_name * ".jl"))
        files = readdir(folder)
        files = map!(file -> joinpath(folder, file), files, files)
        filter!(file -> occursin(module_name * ".jl", file) && occursin(".cov", file), files)
        sources = [read(files[1], String)]
    end
    # Get total coverage for all Julia files
    covered_lines, total_lines = get_summary(coverage)


    clean_folder(LEAP_SRC_PATH)
    clean_folder(TEST_PATH)
    for source in sources
        @info source
    end
    @info "Coverage of $(module_name): $(covered_lines * 100 / total_lines) %"
end


function parse_commandline()
    settings = ArgParseSettings()

    @add_arg_table settings begin
        "--module-coverage", "-c"
            help = "another option with an argument"
            arg_type = String
            default = "all"
    end

    return parse_args(settings)
end


function main()
    parsed_args = parse_commandline()
    process_coverage_files(parsed_args["module-coverage"])
end


main()
