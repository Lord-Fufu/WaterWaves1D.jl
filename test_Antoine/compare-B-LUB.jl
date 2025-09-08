
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

#=
function phi1(x,t)
    return 0.01 * cos.(kc * x) #- wc * t)
end

function phi2(x,t)
    return 0.01 * sin.(kc * x) #- wc * t)
end

modes_list = []=#

phiArray = zeros(Float64, (1, 2))

phiArray[1,1] = kc
phiArray[1,2] = 0.01


#=
dx=L/N_points/2
k_max=1/dx
dk=2/L

#n_phi=N_points
#waveNumbers = LinRange(dk, k_max, n_phi)
n_phi = 1
waveNumbers = [kc]
phiArray = zeros(Float64, (n_phi, 2))

coeff = sum((waveNumbers[1]./ waveNumbers).^5)^(1/2)
total_amplitude = 0.005
max_amplitude = total_amplitude/coeff
amplitudes = max_amplitude * (waveNumbers[1] ./ waveNumbers).^(5/2)

phiArray[:,1] = waveNumbers
phiArray[:,2] = amplitudes
=#

#=for i in 1:n_phi
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
    ϵ  = 0.1,   	            # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = N_points,  	            # number of collocation points
    Ns = 300,
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



z(x) = exp.( - abs.(x).^4);                 # surface deformation
v(x) = zero(x);                             # zero initial velocity
init = Init(z,v);                           # generate the initial data with correct type

B_model=Boussinesq(param)
LUB_model=LUBoussinesq(param)   
WW_model=WaterWaves(param)                  # The water waves system


B_problem=Problem(B_model, init, param ; solver = RK4(B_model))
LUB_problem=LUProblem(LUB_model, init, param ; solver = RK4Sto(LUB_model))
WW_problem=Problem(WW_model, init, param) ;


solve!(B_problem;verbose=true);
LUsolve!(LUB_problem;verbose=true);
convLUB_problem = Problem(LUB_problem)
solve!(WW_problem;verbose=true);




using Plots
plot([B_problem, convLUB_problem, WW_problem], T=5, labels = ["Boussinesq" "LU Boussinesq" "Water Waves"],
xlabel = "Horizontal domain (x)", ylabel = "Surface elevation (η)",
titlefontsize=14,
guidefontsize=14,
tickfontsize=10,
legendfontsize=10,
xlim = (-20,20),
ylim = (-0.4,1.0))

savefig(".//Documents//__plot Julia//plotsForPaper//LUB.pdf")


#energy conservation


B_energyL1_eta = zeros(Float64, param.Ns)
B_energyL2_eta = zeros(Float64, param.Ns)
B_energyL2_velocity = zeros(Float64, param.Ns)

LUB_energyL1_eta = zeros(Float64, param.Ns)
LUB_energyL2_eta = zeros(Float64, param.Ns)
LUB_energyL2_velocity = zeros(Float64, param.Ns)

WW_energyL1_eta = zeros(Float64, param.Ns)
WW_energyL2_eta = zeros(Float64, param.Ns)
WW_energyL2_velocity = zeros(Float64, param.Ns)

using FFTW
for time_index in 1:param.Ns
    B_energyL1_eta[time_index] = sum(real.(ifft(B_problem.data.U[time_index][:,1])))*L/N_points
    B_energyL2_eta[time_index] = sum(abs.(B_problem.data.U[time_index][:,1]).^2)*L/N_points
    B_energyL2_velocity[time_index] = sum(abs.(B_problem.data.U[time_index][:,2]).^2)*L/N_points

    LUB_energyL1_eta[time_index] = sum(real.(ifft(convLUB_problem.data.U[time_index][:,1])))*L/N_points
    LUB_energyL2_eta[time_index] = sum(abs.(convLUB_problem.data.U[time_index][:,1]).^2)*L/N_points
    LUB_energyL2_velocity[time_index] = sum(abs.(convLUB_problem.data.U[time_index][:,2]).^2)*L/N_points

    WW_energyL1_eta[time_index] = sum(real.(ifft(WW_problem.data.U[time_index][:,1])))*L/N_points
    WW_energyL2_eta[time_index] = sum(abs.(WW_problem.data.U[time_index][:,1]).^2)*L/N_points
    WW_energyL2_velocity[time_index] = sum(abs.(WW_problem.data.U[time_index][:,2]).^2)*L/N_points
end

using Plots

ratio = (LUB_energyL2_velocity[300] - LUB_energyL2_velocity[150]) / (LUB_energyL1_eta[300] - LUB_energyL1_eta[150])

timeValues = LinRange(0,param.T,param.Ns)
plot(timeValues, ratio*(B_energyL1_eta), label = "L1 norm eta")
plot!(timeValues, B_energyL2_velocity,  label = "L2 norm u")
plot!(timeValues, B_energyL2_velocity .+ ratio*(B_energyL1_eta),  label = "total energy")

plot(timeValues, ratio*(LUB_energyL1_eta), label = "L1 norm eta")
plot!(timeValues, LUB_energyL2_velocity,  label = "L2 norm u")
plot!(timeValues, LUB_energyL2_velocity .+ ratio*(LUB_energyL1_eta),  label = "total energy")

plot(timeValues, ratio*(WW_energyL1_eta), label = "L1 norm eta")
plot!(timeValues, WW_energyL2_velocity,  label = "L2 norm u")
plot!(timeValues, WW_energyL2_velocity .+ ratio *(WW_energyL1_eta),  label = "total energy")





#n_show = 500
anim = @animate for t in LinRange(0,param.T,param.Ns)
    plot([B_problem, convLUB_problem], T = t)
    ylims!(-0.5, 1)
    xlims!(-50, 50)
end

gif( anim, ".//Documents//__plot Julia//LUB-000.gif", fps = 30 )




plot(smoother, xlims=(-L,L),tickfontsize = 12, legendfontsize = 14, guidefontsize = 14, label = "s_α", xlabel = "Horizontal domain (x)", ylabel = "Smoothing factor")
savefig(".//Documents//__plot Julia//plotsForPaper//smoothingFactor.pdf")
