export RK4Sto
export STOstep!

"""
    RK4(arguments;realdata)

Explicit Runge-Kutta fourth order solver.

Construct an object of type `TimeSolver` to be used in `Problem(model, initial, param; solver::TimeSolver)`

Arguments can be either
0. an object of type `AbstractModel`;
1. an `Array` of size `(N,datasize)` where `N` is the number of collocation points and `datasize` the number of equations solved;
2. `(param,datasize)` where `param is a `NamedTuple` containing a key `N`, and `datasize` a integer (optional, by default `datasize=2`).

The keyword argument `realdata` is optional, and determines whether pre-allocated vectors are real- or complex-valued.
By default, they are either determined by the model or the type of the array in case `0.` and `1.`, complex-valued in case `2.`.

"""
struct RK4Sto <: TimeSolver

    U1 :: Array
    dU :: Array
    label :: String

    function RK4Sto( U :: Array; realdata=nothing )
        U1 = copy(U)
        dU = copy(U)
        if realdata==true
            U1 = real.(U1);dU = real.(dU)
        end
        if realdata==false
            U1 = complex.(U1);dU = complex.(dU)
        end
        new( U1, dU, "stochastic RK4")
    end

    function RK4Sto( model :: AbstractModel; realdata=nothing )
        U=model.mapto(Init(x->0*x,x->0*x))
        RK4Sto( U; realdata=realdata)
    end
    function RK4Sto( param::NamedTuple, datasize=2::Int; realdata=nothing )
        RK4Sto( zeros(Complex{Float64}, (param.N,datasize)) ; realdata=realdata)
    end
end


function STOstep!(s  :: RK4Sto,
                m :: AbstractModel,
                U, dt, currentTime ; realData = true)

    isTrueStrato = (m.isTrueStrato)
    dBt = randn(Complex{Float64}, 2*m.n_phi) .* √(dt)

    s.U1 .= U
    m.f!( s.U1 , currentTime )
    s.dU .= s.U1

    s.U1 .= U .+ dt/2 .* s.U1
    m.f!( s.U1 , currentTime )
    s.dU .+= 2 .* s.U1

    s.U1 .= U .+ dt/2 .* s.U1
    m.f!( s.U1 , currentTime )
    s.dU .+= 2 .* s.U1

    s.U1 .= U .+ dt .* s.U1
    m.f!( s.U1 , currentTime)
    s.dU .+= s.U1

    stoch_incr_input = m.stoch_incr!(U,currentTime)
    
    (N,n_var,_) = size(stoch_incr_input)

    stoch_incr_tot = zeros(Complex{Float64}, (N, n_var))

    for k in 1:(2*m.n_phi)
        stoch_incr_tot[:,1] .+= dBt[k] .* stoch_incr_input[:,1,k]
        stoch_incr_tot[:,2] .+= dBt[k] .* stoch_incr_input[:,2,k]
    end

    if isTrueStrato        
        U_temp = copy(U)
        U_temp .+= stoch_incr_tot

        stoch_incr_input_bis = m.stoch_incr!(U_temp,currentTime+dt)
        stoch_incr_tot_bis = zeros(Complex{Float64}, (N, n_var))

        for k in 1:(2*m.n_phi)
            stoch_incr_tot_bis[:,1] .+= dBt[k] .* stoch_incr_input_bis[:,1,k]
            stoch_incr_tot_bis[:,2] .+= dBt[k] .* stoch_incr_input_bis[:,2,k]
        end

        U .+= dt/6 .* s.dU .+ (stoch_incr_tot .+ stoch_incr_tot_bis)/2
    else
        U .+= dt/6 .* s.dU .+ stoch_incr_tot
    end
end
