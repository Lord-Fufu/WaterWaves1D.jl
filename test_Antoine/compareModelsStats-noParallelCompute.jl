
using WaterWaves1D

wc=2 * π / 10
kc=2 * π / 100
L=50
N_points = 2^11

alpha=10
function smoother(x)                        # Smoothing function (to enforce boundary conditions)
    if abs(x/L) < 1
        return exp(1/alpha^2 * (1-1/(1-(x/L)^2)))
    else
        return 0.0
    end
end

n_begin = 1
n_end = 130

#=n_phi=11
waveNumbers = LinRange(2 * π/20, 2 * π/5, n_phi)
phiArray = zeros(Float64, (n_phi, 2))
n_waves = length(waveNumbers)

coeff = sum((waveNumbers[1]./ waveNumbers).^5)^(1/2)
total_amplitude = 0.1
max_amplitude = total_amplitude/coeff
amplitudes = max_amplitude * (waveNumbers[1] ./ waveNumbers).^(5/2)

phiArray[:,1] = waveNumbers
phiArray[:,2] = amplitudes=#


phiArray = zeros(Float64, (1, 2))

phiArray[1,1] = kc
phiArray[1,2] = 0.005


#=function phi1(x,t)
    return 0.01 * cos.(kc * x) #- wc * t)
end

function phi2(x,t)
    return 0.01 * sin.(kc * x) #- wc * t)
end=#


#=modes_list = []

for i in 1:n_waves
    function temp1(x,t)
        return max_amplitude * (waveNumbers[1]/ waveNumbers[i])^(5/2) * cos.(waveNumbers[i] * x - waveNumbers[i] * t) .* smoother(x)
    end
    function temp2(x,t)
        return max_amplitude * (waveNumbers[1]/ waveNumbers[i])^(5/2) * sin.(waveNumbers[i] * x - waveNumbers[i] * t) .* smoother(x)
    end
    push!(modes_list,temp1)
    push!(modes_list,temp2)
end

function phi_function(x,t)
    N = length(modes_list)
    L=[]
    for i in 1:N
        push!(L, modes_list[i](x,t))
    end
    return L
end=#


param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     	            # shallow-water dimensionless parameter
    ϵ  = 0.01,   	            # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = N_points,  	            # number of collocation points
    Ns = 10,
    L  = L,    	                # half-length of the numerical tank (-L,L)
    T  = 5,     	            # final time of computation
    dt = 0.005,  	            # timestep
    phiArray = phiArray, 
    smoother = smoother,
    #phi = (x,t) -> [phi1(x,t) * smoother(x), phi2(x,t) * smoother(x)],
    isTrueStrato = true,
    noiseType = "Waves",
    additiveNoise = false,
);




z(x) = exp.(-abs.(x).^4);               # surface deformation
v(x) = zero(x);                         # zero initial velocity
init = Init(z,v);               # generate the initial data with correct type

det_model=SerreGreenNaghdi(param) # The water waves system
WW_model=WaterWaves(param)

det_problem=Problem(det_model, init, param; solver = RK4(det_model)) ;
WW_problem=Problem(WW_model, init, param) ;

solve!(det_problem;verbose=false);
solve!(WW_problem;verbose=true);

using DelimitedFiles

print("Solving ", n_end - n_begin + 1 ," LU problems, from ", n_begin, " to ", n_end ,", to compute mean and STD. \n")
print("\n")
for p in n_begin:n_end
    print("Problem ", p, " between ", n_begin, " and ", n_end,". \n")
    LU_model=LUSerreGreenNaghdi(param)
    LU_problem="Anything, just need to erase the previous variable, otherwise simulations are messed up..."
    LU_problem=LUProblem(LU_model, init, param; solver = RK4Sto(LU_model))
    LUsolve!(LU_problem;verbose=true);
    convLU_problem = Problem(LU_problem)
    writedlm(".//Documents//__plot Julia//trajectoriesForStats//__other-SGN-files-p3-smallNoise//traj" * string(p) * ".csv",
                                                                            convLU_problem.data.U[11][:,1])
end
