export LUAiry

"""
	Airy(param;mesh,label)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the linear ([Airy](https://en.wikipedia.org/wiki/Airy_wave_theory)) water waves equations.

# Arguments
- `param::NamedTuple` must contain
    - the shallowness parameter `μ` (set the infinite-layer case if `μ=Inf`);
    - optionally, `ν` the shallow/deep water scaling factor. By default, `ν=1` if `μ≦1` and `ν=1/√μ` otherwise;
    - numerical parameters to construct the mesh of collocation points, if `mesh` is not provided.
- `mesh  :: Mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`.
- `label :: String`: a label for future references (default is `"linear (Airy)"`).


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `Airy.f!` to be called in explicit time-integration solvers;
2. a function `Airy.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `Airy.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
    - `η` is the values of surface deformation at collocation points `x`;
    - `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct LUAiry <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String
	stoch_incr! :: Function
	n_phi   :: Int
	isTrueStrato :: Bool


    function LUAiry(param::NamedTuple; # param is a NamedTuple containing all necessary parameters
		IL	    = false,
		label = "BV-linear (LU-Airy, with nonlinear noise parameter ϵ)"  # using a keyword argument allows the user to supersede the default label.
		 )

		# Set up
		μ = param.μ
		ϵ = param.ϵ
		if !in(:ν,keys(param))
			if μ > 1
				ν = 1/sqrt(μ)
				nu = "1/√μ (deep water case)"
			else
				ν = 1
				nu = "1 (shallow water case)"
			end
		else
			ν = param.ν
			nu = "$ν"
		end
		if μ == Inf || ν==0 || IL == true # infinite layer case
			IL = true;  # IL (=Infinite layer) is a flag to be used thereafter
			μ = 1; ν = 1; # Then we should set μ=ν=1 in subsequent formula.
		end

		if in(:isTrueStrato,keys(param))
			isTrueStrato = param.isTrueStrato
		else
			isTrueStrato = true
		end

		mesh = Mesh(param)
		phi = param.phi
		constPhi = param.constPhi
		n_phi = size(phi(0,0))[1]

		# Print information
		info = "$label model.\n"
		if IL == true
			info *= "└─Infinite depth case.\n"
		else
			info *= "└─Shallowness parameter μ=$μ, \
					scaling parameter ν=$nu.\n"
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."
		
		# Collocation points and Fourier modes
		x, k = mesh.x, mesh.k
		# Fourier multipliers
    	∂ₓ	 = 1im * k            # Differentiation
		if IL == true
			∂ₓF₁ = 1im * sign.(k)
		else
			∂ₓF₁ = 1im * tanh.(√μ*k)/(√μ*ν)
		end
		# Pre-allocation
		fftη = zeros(Complex{Float64}, mesh.N)
        fftv = zeros(Complex{Float64}, mesh.N)
		h = zeros(Complex{Float64}, mesh.N)

		# Evolution equations are ∂t U = f(U)
		function f!(U, currentTime)

			fftη .= U[:,1]
			fftv .= U[:,2]

			h .= 1 .+ ϵ .* real(ifft(fftη))

			tablePhi = zeros(Complex{Float64}, (mesh.N, n_phi))
			stochVisc = zeros(Complex{Float64}, mesh.N)
			us = zeros(Complex{Float64}, mesh.N)

			for index in 1:mesh.N
				tablePhi[index,:] .= phi(x[index],currentTime)[:]
			end
			for k in 1:n_phi
				stochVisc .+= tablePhi[:,k] .^ 2
			end

			us .= ϵ .^2 ./ 2 .* ifft(∂ₓ .* fft(stochVisc))
			stochVisc .+= constPhi(currentTime) .^ 2

		   	U[:,1] .= -∂ₓF₁.*fftv 	.+ ϵ .* ∂ₓ .* fft(h .* us )
		   	U[:,2] .= -∂ₓ.*fftη 	.+ ϵ .* fft(us .* ifft(∂ₓ .* fftv))
			if !(isTrueStrato)
				U[:,1] .+= ϵ .^ 2 ./ 2 .* ∂ₓ .* fft(stochVisc .* ifft(∂ₓ .* fftη))
				U[:,2] .+= ϵ .^ 2 ./ 2 .* ∂ₓ .* fft(stochVisc .* ifft(∂ₓ .* fftv))
			end

		end

		# Build raw data from physical data (discrete Fourier transform)
		function mapto(data::InitialData)
			U = [fft(data.η(x)) fft(data.v(x))]
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2])),mesh.x
		end

		function stoch_incr!(U,currentTime)
			
			fftη .= U[:,1]
			fftv .= U[:,2]

			tablePhi = zeros(Complex{Float64}, (mesh.N, n_phi))
			stoch_incr_U = zeros(Complex{Float64},(mesh.N, 2, n_phi+1))

			for i in 1:mesh.N
				tablePhi[i,:] .= phi(x[i],currentTime)[:]
			end

			stoch_incr_U[:,1,1] .-= ϵ .* constPhi.(currentTime) .* ∂ₓ .* fftη[:]
			stoch_incr_U[:,2,1] .-= ϵ .* constPhi.(currentTime) .* ∂ₓ .* fftv[:]
			for k in 1:n_phi
				stoch_incr_U[:,1,k+1] .-= ϵ .* fft( tablePhi[:,k] .* ifft(∂ₓ .* fftη))[:]
				stoch_incr_U[:,2,k+1] .-= ϵ .* fft( tablePhi[:,k] .* ifft(∂ₓ .* fftv))[:]
			end

			return stoch_incr_U
		end

		new(label, f!, mapto, mapfro, info, stoch_incr!, n_phi, isTrueStrato )

    end
end
