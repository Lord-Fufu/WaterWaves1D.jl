export LUEulerMaruyama
export STOstep!

"""
    Euler(arguments;realdata)

Explicit Euler-Maruyama solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,datasize)` where `N` is the number of collocation points and `datasize` the number of equations solved;
2. `(param,datasize)` where `param is a `NamedTuple` containing a key `N`, and `datasize` a integer (optional, by default `datasize=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued in case `2.`.

"""
struct LUEulerMaruyama <: TimeSolver

    U1 :: Array
    label :: String

    function LUEulerMaruyama( U :: Array ; realdata=nothing )
        U1 = copy(U)
        if realdata==true
            U1 = real.(U1);
        end
        if realdata==false
            U1 = complex.(U1);
        end
        new( U1, "Euler-Maruyama" )
    end

    function LUEulerMaruyama( model :: AbstractModel; realdata=nothing )
        U=model.mapto(Init(x->0*x,x->0*x))
        LUEulerMaruyama( U ; realdata=realdata)
    end
    function LUEulerMaruyama( param::NamedTuple, datasize=2::Int; realdata=nothing )
        LUEulerMaruyama( zeros(Complex{Float64}, (param.N,datasize)); realdata=realdata)
    end
end

function STOstep!(solver :: LUEulerMaruyama,
                model :: AbstractModel,
                U , dt, currentTime ; realData = true)
    
    isTrueStrato = (model.isTrueStrato)
    dBt = randn(model.n_phi) .* √(dt)

    solver.U1 .= U
    model.f!( solver.U1, currentTime ) 
    
    if realData
        stoch_incr_input = real.(model.stoch_incr!(U,currentTime))
    else
        stoch_incr_input = complex.(model.stoch_incr!(U,currentTime))
    end
    
    (N,n_var,n_phi) = size(stoch_incr_input)

    stoch_incr_tot = zeros(Complex{Float64}, (N, n_var))

    for k in 1:(model.n_phi)
        stoch_incr_tot[:,1] .+= dBt[k] * stoch_incr_input[:,1,k]
        stoch_incr_tot[:,2] .+= dBt[k] * stoch_incr_input[:,2,k]
    end

    if isTrueStrato        
        U_temp = copy(U)
        U_temp .+= stoch_incr_tot

        if realData
            stoch_incr_input_bis = real.(model.stoch_incr!(U_temp,currentTime+dt))
        else
            stoch_incr_input_bis = complex.(model.stoch_incr!(U_temp,currentTime+dt))
        end
        stoch_incr_tot_bis = zeros(Complex{Float64}, (N, n_var))

        for k in 1:(model.n_phi)
            stoch_incr_tot_bis[:,1] .+= dBt[k] * stoch_incr_input_bis[:,1,k]
            stoch_incr_tot_bis[:,2] .+= dBt[k] * stoch_incr_input_bis[:,2,k]
        end
        U .+= dt .* solver.U1 + (stoch_incr_tot + stoch_incr_tot_bis)/2
    else
        U .+= dt .* solver.U1 + stoch_incr_tot
    end

end
