abstract type Simulation_Module end

abstract type Agent_Module <: Simulation_Module end
abstract type Demographics_Module <: Simulation_Module end
abstract type AsthmaOccurrence_Module <: Simulation_Module end
abstract type AsthmaOutcomes_Module <: Simulation_Module end
abstract type RiskFactors_Module <: Simulation_Module end
abstract type Payoffs_Module <: Simulation_Module end

# demographics
abstract type Birth_Module <: Demographics_Module  end
abstract type Immigration_Module <: Demographics_Module  end
abstract type Emigration_Module <: Demographics_Module  end
abstract type Death_Module <: Demographics_Module  end

# risk factors
abstract type AntibioticExposure_Module   <: RiskFactors_Module  end
abstract type FamilyHistory_Module <: RiskFactors_Module end

# asthma occurrence
abstract type Incidence_Module   <: AsthmaOccurrence_Module  end
abstract type Diagnosis_Module <: AsthmaOccurrence_Module end
abstract type Reassessment_Module <: AsthmaOccurrence_Module end

# asthma outcomes
abstract type Exacerbation_Module   <: AsthmaOutcomes_Module  end
abstract type ExacerbationHist_Module end
abstract type Exacerbation_Severity_Module <: Exacerbation_Module end
abstract type Control_Module   <: AsthmaOutcomes_Module  end

# payoffs
abstract type Utility_Module <: Payoffs_Module end
abstract type Cost_Module <: Payoffs_Module end
