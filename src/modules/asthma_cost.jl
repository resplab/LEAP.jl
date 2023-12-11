struct Cost <: Cost_Module
    parameters
end

function process(ag::Agent,cost::Cost)
    if !ag.has_asthma
        # no asthma, so return 0
        return 0
    else
        # costs due to having exacerbation and control
        return sum(ag.exac_sev_hist.current_year .* cost.parameters[:exac]) + sum(ag.control .* cost.parameters[:control])
    end
end
