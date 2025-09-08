
"""Parameters for (LU)Boussinesq and WaterWaves models. The noise sigma.dB is decomposed as cst.dB0 + phi1.dB1 + phi2.dB2,
where phi1 and phi2 are progessive waves.  """

using WaterWaves1D

wc=2 * π / 10       # Angular frequency of sigma
kc=2 * π / 10       # Wave vector of sigma
L=50               # Half length of the numrical tank
N_points=2^11

alpha=10
function smoother(x)                        # Smoothing function (to enforce boundary conditions)
    if abs(x/L) < 1
        return exp(1/alpha^2 * (1-1/(1-(x/L)^2)))
    else
        return 0.0
    end
end


dx=L/N_points/2
k_max=1/dx
dk=2/L

n_phi=N_points
waveNumbers = LinRange(dk, k_max, n_phi)
phiArray = zeros(Float64, (n_phi, 2))

coeff = sum((waveNumbers[1]./ waveNumbers).^5)^(1/2)
total_amplitude = 0.005
max_amplitude = total_amplitude/coeff
amplitudes = max_amplitude * (waveNumbers[1] ./ waveNumbers).^(5/2)

phiArray[:,1] = waveNumbers
phiArray[:,2] = amplitudes


param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     	            # shallow-water dimensionless parameter
    ϵ  = 0.1,   	            # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = N_points,  	            # number of collocation points
    Ns = 150,
    L  = L,    	                # half-length of the numerical tank (-L,L)
    T  = 5,     	            # final time of computation
    dt = 0.005,  	            # timestep
    phiArray = phiArray, 
    smoother = smoother,
    #phi = (x,t) -> [phi1(x,t) * smoother(x), phi2(x,t) * smoother(x)],
    isTrueStrato = false,
    noiseType = "Waves",
);

paramWW = param #merge(param,(L=20, dt=0.01, T=15))        #params for WaterWaves model

z(x) = exp.( - abs.(x).^4);                 # surface deformation
v(x) = zero(x);                             # zero initial velocity
init = Init(z,v);                           # generate the initial data with correct type
