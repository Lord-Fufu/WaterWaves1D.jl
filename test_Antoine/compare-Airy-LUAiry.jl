

using WaterWaves1D

wc=2 * π / 10
kc=2 * π / 25
L=100

function smoother(x)
    if abs(x/L) < 1
        return exp(-1/(1-(x/L)^2) )
    else
        return 0
    end
end

function constPhi(t)
    return 0 #0.01
end

function phi1(x,t)
    return 0.05 * cos.(kc * x) #- wc * t)
end

function phi2(x,t)
    return 0.05 * sin.(kc * x) #- wc * t)
end

param = (
    # Physical parameters. Variables are non-dimensionalized as in Lannes, The water waves problem, isbn:978-0-8218-9470-5
    μ  = 1,     	            # shallow-water dimensionless parameter
    ϵ  = 3/4,   	            # nonlinearity dimensionless parameter
    # Numerical parameters
    N  = 2^11,  	            # number of collocation points
    L  = L,    	            # half-length of the numerical tank (-L,L)
    T  = 20,     	            # final time of computation
    dt = 0.01,  	            # timestep
    constPhi = t -> constPhi(t),
    phi= (x,t) -> [phi1(x,t) * smoother(x), phi2(x,t) * smoother(x)],            # noise data
    isTrueStrato = false
                );


z(x) = exp.( - abs.(x).^4);               # surface deformation
v(x) = zero(x);                         # zero initial velocity
init = Init(z,v);               # generate the initial data with correct type

Airy_model=Airy(param)
LUAiry_model=LUAiry(param)   # The water waves system

Airy_problem=Problem(Airy_model, init, param ; solver = Euler(Airy_model))
LUAiry_problem=LUProblem(LUAiry_model, init, param)

solve!(Airy_problem;verbose=true);
LUsolve!(LUAiry_problem;verbose=true);
convLUAiry_problem = Problem(LUAiry_problem)



using Plots
plot([Airy_problem,convLUAiry_problem])

n_show = 500
anim = @animate for t in LinRange(0,param.T,n_show)
    plot([Airy_problem,convLUAiry_problem], T = t)
    ylims!(-0.5, 1)
    xlims!(-20,20)
end

gif(anim, ".//Documents//__plot Julia//LU-Airy.gif", fps=round(Int, n_show ./ param.T))


