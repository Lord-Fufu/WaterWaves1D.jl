
using WaterWaves1D

wc=2 * π / 10
kc=2 * π / 100
L=50
L_smooth = 50
N_points = 2^11


alpha_smooth=10
function smoother(x)                        # Smoothing function (to enforce boundary conditions)
    if abs(x/L_smooth) < 1
        return exp(1/alpha_smooth^2 * (1-1/(1-(x/L_smooth)^2)))
    else
        return 0.0
    end
end


#=function phi1(x,t)
    return 0.001 * cos.(kc * x) #- wc * t)
end

function phi2(x,t)
    return 0.001 * sin.(kc * x) #- wc * t)
end=#

phiArray = zeros(Float64, (1, 2))

phiArray[1,1] = kc
phiArray[1,2] = 0.01

#=
dx=L/N_points/2
k_max=1/dx
dk=2/L

n_phi=N_points
waveNumbers = LinRange(dk, k_max, n_phi)
#n_phi = 1
#waveNumbers = [kc]
phiArray = zeros(Float64, (n_phi, 2))

coeff = sum((waveNumbers[1]./ waveNumbers).^5)^(1/2)
total_amplitude = 0.001
max_amplitude = total_amplitude/coeff
amplitudes = max_amplitude * (waveNumbers[1] ./ waveNumbers).^(5/2)

phiArray[:,1] = waveNumbers
phiArray[:,2] = amplitudes
=#

param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 0.01,     	            # shallow-water dimensionless parameter
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


z(x) = exp.( - abs.(x).^4);               # surface deformation
#=L_init = 2
alpha_init = 1
function initialProfile(x)
    if abs(x/L_init) < 1.0
        return exp(1/alpha_init^2 * (1-1/(1-(x/L_init)^2)))
    else
        return 0.0
    end
end
z(x) = initialProfile.(x)=#
v(x) = zero(x);                         # zero initial velocity
init = Init(z,v);               # generate the initial data with correct type



#=using Plots
plot(smoother, -L, L, label = "s_α", xlabel = "Horizontal domain (x)", ylabel = "Smoothing factor",
titlefontsize=14,
guidefontsize=14,
tickfontsize=10,
legendfontsize=10,
xlim = (-L-5,L+5),
ylim = (-0.1,1.1))
savefig(".//Documents//__plot Julia//plotsForPaper//BCfactor.pdf")

plot(z, -L, L, label = "Initial surface elevation", xlabel = "Horizontal domain (x)", ylabel = "Surface elevation (η)",
titlefontsize=14,
guidefontsize=14,
tickfontsize=10,
legendfontsize=10,
xlim = (-L-5,L+5),
ylim = (-0.1,1.1))
savefig(".//Documents//__plot Julia//plotsForPaper//initialCondition.pdf")=#





SV_model=SaintVenant(param)
LUSV_model=LUSaintVenant(param)   # The water waves system
WW_model=WaterWaves(param)


SV_problem=Problem(SV_model, init, param ; solver = RK4(SV_model))
LUSV_problem=LUProblem(LUSV_model, init, param ; solver = RK4Sto(LUSV_model))
WW_problem=Problem(WW_model, init, param);


solve!(SV_problem;verbose=true);
LUsolve!(LUSV_problem;verbose=true);
convLUSV_problem = Problem(LUSV_problem);
solve!(WW_problem;verbose=true);


using Plots
plot([SV_problem, convLUSV_problem, WW_problem], T=5, labels = ["Saint-Venant" "LU Saint-Venant" "Water Waves"],
xlabel = "Horizontal domain (x)", ylabel = "Surface elevation (η)",
titlefontsize=14,
guidefontsize=14,
tickfontsize=10,
legendfontsize=10,
xlim = (-20,20),
ylim = (-0.4,1.0))




savefig(".//Documents//__plot Julia//plotsForPaper//LUSV.pdf")



#-------Energy conservation-------#

SV_m  = zeros(Float64, param.Ns)
SV_mom = zeros(Float64, param.Ns)
SV_energ = zeros(Float64, param.Ns)

LUSV_m  = zeros(Float64, param.Ns)
LUSV_mom = zeros(Float64, param.Ns)
LUSV_energ = zeros(Float64, param.Ns)

#WW_m  = zeros(Float64, param.Ns)
#WW_mom = zeros(Float64, param.Ns)
#WW_energ = zeros(Float64, param.Ns)

function momentumdiff_test(p::Problem, ϵ::Float64; T=nothing,rel=false)
	η,mom,x = solution(p;T=T)
	η0,mom0,x0 = solution(p;T=0)
	if !(x[2:end].-x[2]≈x[1:end-1].-x[1])
		@error("The horizontal impulse difference cannot be computed because the solution is defined on a non-regularly spaced mesh.")
	else
        if rel==false return sum(mom - mom0)*(x[2]-x[1]) else return sum((1 .+ ϵ*η).*v .- (1 .+ ϵ*η0).*v0)/sum((1 .+ ϵ*η0).*v0) end
	end
end

for time_index in 1:param.Ns
    t = time_index*param.T/param.Ns
    #WW_m[time_index]   = massdiff(WW_problem, T=t)
    SV_m[time_index]    = massdiff(SV_problem, T=t)
    LUSV_m[time_index]  = massdiff(convLUSV_problem, T=t)

    #WW_mom[time_index]  = momentumdiff(WW_problem, T=t)
    SV_mom[time_index]   = momentumdiff(SV_problem, T=t)
    LUSV_mom[time_index] = momentumdiff_test(convLUSV_problem, param.ϵ, T=t)

    #WW_energ[time_index]   = energydiff(WW_problem, T=t)
    #SV_energ[time_index]   = energydiff(SV_problem, T=t)
    #LUSV_energ[time_index] = energydiff(convLUSV_problem, T=t)
end

using Plots

timeValues = LinRange(0,param.T,param.Ns)

plot(timeValues, SV_m, label = "Saint-Venant Mass")
plot!(timeValues, LUSV_m, label = "LU Saint-Venant Mass")

plot(timeValues, SV_mom, label = "Saint-Venant Momentum")
plot!(timeValues, LUSV_mom, label = "LU Saint-Venant Momentum")

plot(timeValues, SV_energ, label = "Saint-Venant Energy")
plot!(timeValues, LUSV_energ, label = "LU Saint-Venant Energy")

n_show = 300
anim = @animate for t in LinRange(0,param.T,n_show)
    plot([SV_problem, convLUSV_problem], T = t)
    ylims!(-0.5, 1)
    xlims!(-50,50)
end

gif(anim, "C://Users//amoneyro//Documents//__plot Julia//LUSV-test.gif", fps=round(Int, n_show ./ param.T))



n_show = 300
anim = @animate for t in LinRange(0,param.T,n_show)
    eta,vel,space = solution(convLUSV_problem, T=t)
    plot(space, eta .* vel)
    plot!(space, -reverse(eta .* vel))
    ylims!(-0.5, 0.5)
    xlims!(-L,L)
end

gif(anim, "C://Users//amoneyro//Documents//__plot Julia//LUSV-mom.gif", fps=round(Int, n_show ./ param.T))







using FFTW

mesh=Mesh(param)
∂ₓ = mesh.k

phiArray = param.phiArray
smoother = param.smoother
smoothingField = smoother.(mesh.x)
smoothingField2 = smoothingField.^2
stochVisc = sum(phiArray[:,2].^2) * smoothingField2
fftus = param.ϵ .^2 ./ 2 .* ∂ₓ .* fft(stochVisc)
us = realRounder.(real.(ifft(fftus)))

plot(us, -L, L)

N=N_points
temp_fftSigma = zeros(N)
temp_fftSigma[1+round(Int, phiArray[1,1]/dk)] = phiArray[1,2]

temp_sigma1 = N*(real.(ifft(temp_fftSigma)) .* smoothingField .+ (1 .- smoothingField) .* reverse(real.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )
temp_sigma2 = N*(smoothingField .* imag.(ifft(temp_fftSigma))) #.+ (1 .- smoothingField) .* reverse(imag.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )

