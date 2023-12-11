struct Utility <: Utility_Module
    parameters
end

function process(ag::Agent,util::Utility)
    # baseline util
    baseline = util.parameters[:eq5d][(ag.age,Int(ag.sex))].eq5d[1]
    if !ag.has_asthma
        # no asthma, so return the baseline value
        return baseline
    else
        # disutility due to having exacerbation and control
        return max(0,baseline - sum(ag.exac_sev_hist.current_year .* util.parameters[:exac]) - sum(ag.control .* util.parameters[:control]))
    end
end
