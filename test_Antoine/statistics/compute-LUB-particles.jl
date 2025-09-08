
@everywhere path = "C:\\Users\\amoneyro\\.julia\\packages\\WaterWaves1D\\1xDqi\\test_Antoine\\"
@everywhere include(path * "statistics\\compute-LUB-param.jl")
#@everywhere include(path * "statistics\\compute-LUB-param-freqBand.jl")


#Define the deterministic models (on the master worker, i.e worker 1)
#Define the stochastic model (on every worker)
B_model=Boussinesq(param)
WW_model=WaterWaves(paramWW)

@everywhere using DistributedArrays



#Define the problems on solve them (on the master worker)
B_problem=Problem(B_model, init, param ; solver = RK4(B_model))
WW_problem=Problem(WW_model, init, paramWW) ;
solve!(B_problem;verbose=true);
solve!(WW_problem;verbose=true);




n_loop = n_sample÷n_processors
n_rest = n_sample%n_processors

LUB_model=LUBoussinesq(param)
#Broadcast information
info_problem=LUProblem(LUB_model, init, param)
@info "Now solving $n_sample problems of type $(info_problem.label) to compute the mean and std,\n\
    with timestep dt=$(info_problem.times.dt), final time T=$(info_problem.times.tfin),\n\
    and N=$(param.N) collocation points."

#Create distributed array for storing the results
loc_array = []
for k in 1:n_sample
    push!(loc_array, Problem(info_problem))
end
pbmList_glob = distribute(loc_array)

#Solve many stochastic problems, and store the results in the distributed array
for p in 0:(n_loop-1)
    @sync begin @distributed for k in 0:(n_processors-1)
        print("Problem ", 1 + k + p*n_processors, " over ", n_sample, ".\n")
        LUB_model=LUBoussinesq(param)
        temp_problem="Anything, just need to erase the previous variable, otherwise simulations are messed up..."
        temp_problem=LUProblem(LUB_model, init, param ; solver = RK4Sto(LUB_model))
        LUsolve!(temp_problem;verbose=false);
        remotecall_wait(D->localpart(D)[1+p] = Problem(temp_problem), myid(), pbmList_glob)
    end end
    print("\n")
end

@sync begin @distributed for k in 0:(n_rest-1)
    print("Problem ", 1+k+n_loop*n_processors, " over ", n_sample, ".\n")
    LUB_model=LUBoussinesq(param)
    temp_problem="Anything, just need to erase the previous variable, otherwise simulations are messed up..."
    temp_problem=LUProblem(LUB_model, init, param ; solver = RK4Sto(LUB_model))
    LUsolve!(temp_problem;verbose=false);
    remotecall_wait(D->localpart(D)[1+n_loop] = Problem(temp_problem), myid(), pbmList_glob)
end end

print("Solving done. Now converting distributed array into a non-distributed one. \n")

#Convert the distributed array into a vector of which elements have the type Problem
pbmList_glob = convert(Vector{Any}, pbmList_glob)

print("Done.")