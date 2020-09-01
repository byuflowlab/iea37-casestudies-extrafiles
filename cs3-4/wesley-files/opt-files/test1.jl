layout_number = Base.parse(Int, ENV["SLURM_ARRAY_TASK_ID"])
println("First argument is: ", Base.ARGS[1])
println("Second argument is: ", Base.ARGS[2])
println("Layout number is: ", layout_number)