abstract type SimulationModule end

abstract type AgentModule <: SimulationModule end
abstract type DemographicsModule <: SimulationModule end
abstract type AsthmaOccurrenceModule <: SimulationModule end
abstract type AsthmaOutcomesModule <: SimulationModule end
abstract type RiskFactorsModule <: SimulationModule end
abstract type PayoffsModule <: SimulationModule end

# demographics
abstract type BirthModule <: DemographicsModule  end
abstract type ImmigrationModule <: DemographicsModule  end
abstract type EmigrationModule <: DemographicsModule  end
abstract type DeathModule <: DemographicsModule  end

# risk factors
abstract type AntibioticExposureModule   <: RiskFactorsModule  end
abstract type FamilyHistoryModule <: RiskFactorsModule end

# asthma occurrence
abstract type IncidenceModule   <: AsthmaOccurrenceModule  end
abstract type DiagnosisModule <: AsthmaOccurrenceModule end
abstract type ReassessmentModule <: AsthmaOccurrenceModule end

# asthma outcomes
abstract type ExacerbationModule   <: AsthmaOutcomesModule  end
abstract type ExacerbationHistModule end
abstract type ExacerbationSeverityModule <: ExacerbationModule end
abstract type ExacerbationSeverityHistModule end
abstract type ControlModule   <: AsthmaOutcomesModule  end

# payoffs
abstract type UtilityModule <: PayoffsModule end
abstract type CostModule <: PayoffsModule end


# pollution
abstract type CensusDivisionModule <: SimulationModule end
abstract type CensusTableModule <: SimulationModule end
abstract type CensusBoundariesModule <: SimulationModule end
abstract type PollutionModule <: SimulationModule end
abstract type PollutionTableModule <: SimulationModule end

abstract type OutcomeMatrixModule <: SimulationModule end
