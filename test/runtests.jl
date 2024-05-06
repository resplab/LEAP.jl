using Test
using LEAP


@testset "All Tests" begin

    if isempty(ARGS) ||  "all" in ARGS
        all_tests = true
    else
        all_tests = false
    end

    if all_tests || "agent" in ARGS
        include("test_agent.jl")
    end
    if all_tests || "antibiotic" in ARGS
        include("test_antibioticexposure.jl")
    end
    if all_tests || "cost" in ARGS
        include("test_asthmacost.jl")
    end
    if all_tests || "birth" in ARGS
        include("test_birth.jl")
    end
    if all_tests || "census" in ARGS
        include("test_censusdivision.jl")
    end
    if all_tests || "control" in ARGS
        include("test_control.jl")
    end
    if all_tests || "death" in ARGS
        include("test_death.jl")
    end
    if all_tests || "emigration" in ARGS
        include("test_emigration.jl")
    end
    if all_tests || "exacerbation" in ARGS
        include("test_exacerbation.jl")
    end
    if all_tests || "family_history" in ARGS
        include("test_familyhistory.jl")
    end
    if all_tests || "immigration" in ARGS
        include("test_immigration.jl")
    end
    if all_tests || "incidence" in ARGS
        include("test_incidence.jl")
    end
    if all_tests || "outcome_matrix" in ARGS
        include("test_outcomematrix.jl")
    end
    if all_tests || "pollution" in ARGS
        include("test_pollution.jl")
    end
    if all_tests || "prevalence" in ARGS
        include("test_prevalence.jl")
    end
    if all_tests || "reassessment" in ARGS
        include("test_reassessment.jl")
    end
    if all_tests || "severity" in ARGS
        include("test_severity.jl")
    end
    if all_tests || "simulation" in ARGS
        include("test_simulation.jl")
    end
    if all_tests || "utility" in ARGS
        include("test_utility.jl")
    end
end
