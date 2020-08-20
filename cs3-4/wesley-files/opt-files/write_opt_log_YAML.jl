"""
    write_opt_log_YAML(file_name, data)

write optimization log and related information to .yaml

# Arguments
- `file_name::String`: path/and/name/of/location/file.yaml
- `funcalls::Array{Array{Float,1}}`: log of all AEP function calls
- `baseyaml::String`: path/and/name/of/location/file.yaml
"""

function write_opt_log_YAML(filename, funcalls; baseyaml=string(@__DIR__, "/default_cs4_log.yaml"), title="", titledescription="", 
    gradient_based=true, algorithm_name="", program_language="Julia", total_optimizations=[], total_wall_time=[], 
    units="s", aepunits="MWh")

    ### Retrieve function call history and auxiliary file names from <.yaml> file.

    # read in the base .yaml file
    base = YAML.load(open(baseyaml))

    # get number of optimizations
    if typeof(funcalls)==Array{Array{Float64,1},1}
        noptimizations = length(funcalls)
    elseif typeof(funcalls)==Array{Float64,1}
        noptimizations = 1
    else
        error("invalid object type for \"funcalls\"")
    end
    println(noptimizations)

    # save the title and description to the yaml database
    base["title"] = title
    base["description"] = titledescription

    # save optimization summary information
    base["optimization_summary"]["algorithm_name"] = algorithm_name
    base["optimization_summary"]["total_wall_time"] = total_wall_time
    base["optimization_summary"]["program_language"] = program_language

    # save function call history
    base["optimization_summary"]["optimization_log_1"]["function_calls"] = length(funcalls[1])
    base["optimization_summary"]["optimization_log_2"]["function_calls"] = length(funcalls[2])
    base["optimization_summary"]["optimization_log_3"]["function_calls"] = length(funcalls[3])

    # save function call history
    base["optimization_summary"]["optimization_log_1"]["annual_energy_production"] = funcalls[1]
    base["optimization_summary"]["optimization_log_2"]["annual_energy_production"] = funcalls[2]
    base["optimization_summary"]["optimization_log_3"]["annual_energy_production"] = funcalls[3]

    # write results to YAML file
    YAML.write_file(filename, base)
end