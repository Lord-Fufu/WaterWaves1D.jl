
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




z(x) = exp.(-abs.(x).^4);               # surface deformation
v(x) = zero(x);                         # zero initial velocity
init = Init(z,v);               # generate the initial data with correct type

SGN_model=SerreGreenNaghdi(param) # The water waves system
LUSGN_model=LUSerreGreenNaghdi(param)   # The water waves system
WW_model=WaterWaves(param)

SGN_problem=Problem(SGN_model, init, param; solver = RK4(SGN_model)) ;
LUSGN_problem=LUProblem(LUSGN_model, init, param; solver = RK4Sto(LUSGN_model))
WW_problem=Problem(WW_model, init, param) ;

solve!(SGN_problem;verbose=false);
LUsolve!(LUSGN_problem;verbose=true);
convLUSGN_problem = Problem(LUSGN_problem)
solve!(WW_problem;verbose=true);


SGN_problem.label="Serre-Green-Naghdi"
convLUSGN_problem.label="LU-SGN"


using Plots
plot([SGN_problem, convLUSGN_problem, WW_problem], T=5, labels = ["SGN" "LU SGN" "Water Waves"],
xlabel = "Horizontal domain (x)", ylabel = "Surface elevation (η)",
titlefontsize=14,
guidefontsize=14,
tickfontsize=10,
legendfontsize=10,
xlim = (-20,20),
ylim = (-0.4,1.0))

savefig(".//Documents//__plot Julia//plotsForPaper//LUSGN.pdf")


n_show = 100
anim = @animate for t in LinRange(0,param.T,n_show)
    plot([SGN_problem,convLUSGN_problem], T = t)
    ylims!(-0.5, 1)
    xlims!(-50,50)
end

gif(anim, ".//Documents//__plot Julia//comparison//LUSGN-Fourier-explode.gif", fps=round(Int, n_show ./ param.T))

