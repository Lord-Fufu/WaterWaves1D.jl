
"""Parameters for (LU)Boussinesq and WaterWaves models. The noise sigma.dB is decomposed as cst.dB0 + phi1.dB1 + phi2.dB2,
where phi1 and phi2 are progessive waves.  """

using WaterWaves1D

wc=2 * π / 10       # Angular frequency of sigma
kc=2 * π / 100       # Wave vector of sigma
L=50               # Half length of the numrical tank
N_points = 2^11

alpha=10
function smoother(x)                        # Smoothing function (to enforce boundary conditions)
    if abs(x/L) < 1
        return exp(1/alpha^2 * (1-1/(1-(x/L)^2)))
    else
        return 0.0
    end
end


#=kc = exp.(LinRange(log(2 * π/10), log(2 * π/0.1), 20) )
wc = kc 

coeff = sum((kc[1]./ kc).^5)^(1/2)
total_amplitude = 0.1
max_amplitude = total_amplitude/coeff

n_waves = length(wc)

n_phi=11
waveNumbers = LinRange(2 * π/20, 2 * π/5, n_phi)
phiArray = zeros(Float64, (n_phi, 2))
n_waves = length(waveNumbers)

coeff = sum((waveNumbers[1]./ waveNumbers).^5)^(1/2)
total_amplitude = 1
max_amplitude = total_amplitude/coeff
amplitudes = max_amplitude * (waveNumbers[1] ./ waveNumbers).^(5/2)

phiArray[:,1] = waveNumbers
phiArray[:,2] = amplitudes

modes_list = []

for i in 1:n_waves
    function temp1(x,t)
        return max_amplitude * (kc[1]/ kc[i])^(5/2) * cos.(kc[i] * x - wc[i] * t) .* smoother(x)
    end
    function temp2(x,t)
        return max_amplitude * (kc[1]/ kc[i])^(5/2) * sin.(kc[i] * x - wc[i] * t) .* smoother(x)
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
end

=#


#=
function phi1(x,t)                          # Oscillating noises phi1 and phi2 (progressive waves)
    return 0.01 * cos.(kc * x) #- wc * t)
end

function phi2(x,t)
    return 0.01 * sin.(kc * x) #- wc * t)
end=#

phiArray = zeros(Float64, (1, 2))

phiArray[1,1] = kc
phiArray[1,2] = 0.001

param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 0.01,     	            # shallow-water dimensionless parameter
    ϵ  = 0.1,   	            # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = N_points,  	            # number of collocation points
    Ns = 5,
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

paramWW = param #merge(param,(L=20, dt=0.01, T=15))        #params for WaterWaves model

z(x) = exp.( - abs.(x).^4);                 # surface deformation
v(x) = zero(x);                             # zero initial velocity
init = Init(z,v);                           # generate the initial data with correct type
                      # generate the initial data with correct type
