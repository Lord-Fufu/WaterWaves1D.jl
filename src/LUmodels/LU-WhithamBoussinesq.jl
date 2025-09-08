export LUWhithamBoussinesq

"""
    WhithamBoussinesq(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
a Boussinesq-type model with full-dispersion property.

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `Boussinesq`: if `true` (default is `false`), compute the standard Boussinesq system instead (see `Boussinesq(param;kwargs)`);
- a parameter `α` which determines the model solved:
    - If `α = 1` (default), then the model has been introduced in [Dinvay, Dutykh and Kalisch](https://doi.org/10.1016/j.apnum.2018.09.016);
    - If `α = 1/2`, then the model is a quasilinear version;
    - If `α < 1/2`, then expect instabilities stemming from ill-posedness of the model.
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `ktol`: tolerance of the low-pass Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Whitham-Boussinesq"`);


# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WhithamBoussinesq.f!` to be called in explicit time-integration solvers;
2. a function `WhithamBoussinesq.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed.
3. a function `WhithamBoussinesq.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
	- `η` is the values of surface deformation at collocation points `x`;
	- `v` is the derivative of the trace of the velocity potential at `x`.

"""
mutable struct LUWhithamBoussinesq <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	info 	:: String
	stoch_incr! :: Function
	n_phi   :: Int
	isTrueStrato :: Bool
	smoothingField :: Vector{Float64}



    function LUWhithamBoussinesq(param::NamedTuple; Boussinesq=false, SGN=false,
								α = 1, a = -1//3, b = 1//3,
								dealias = 0,
								ktol	= 0,
								label 	= nothing,
								)

		# Set up
		μ 	= param.μ
		ϵ 	= param.ϵ
		dt  = param.dt
		additiveNoise = Float64(param.additiveNoise)
		
		if in(:noiseType,keys(param))
			noiseType = param.noiseType
		else
			noiseType = "notWaves"
		end

		if in(:isTrueStrato,keys(param))
			isTrueStrato = param.isTrueStrato
		else
			isTrueStrato = true
		end

		mesh = Mesh(param)


		if Boussinesq == true && SGN == false
			if isnothing(label) label = "LU-Boussinesq" end
			info_param = "a=$a, b=$b, c=0 and d=$b"
		elseif Boussinesq == true && SGN == true
			if isnothing(label) label = "LU-Serre-Green-Naghdi" end
			info_param = "a=$a, b=$b, c=0 and d=$b"
		else
			if isnothing(label) label = "LU-Saint-Venant" end
			info_param = "a=0, b=0, c=0 and d=0"
		end

		if in(:phiArray,keys(param)) && in(:smoother,keys(param))
			phiArray = param.phiArray
			smoother = param.smoother
			smoothingField = smoother.(mesh.x)
			smoothingField2 = realRounder.(smoothingField.^2)
			n_phi = length(phiArray[:,1])
			if noiseType != "Waves"
				@error "The given data is suitable for noise type 'Waves'. Please change the parameter 'noiseType'."
			end
		elseif in(:phi,keys(param))
			phi = param.phi
			n_phi = size(phi(0,0))[1]
			if noiseType == "Waves"
				@error "Please input Fourier wave vector and amplitude with parameter 'phiArray' when using noise type 'Waves'."
			end
		end

		# Print information
		info = "$label model with $info_param.\n"
		info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"
		if dealias == 0
			info *= "└─No dealiasing. "
		else
			info *= "└─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end

		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate data
		k = mesh.k
		x 	= mesh.x
		∂ₓ	=  1im * k

		dk = mesh.dk
		N  = mesh.N

		if Boussinesq==false
			F₁ 	= 1
			F₂ = 1
		else
			F₂ = 1 ./(1 .+ μ*b*abs.(k).^2)
			F₁ 	= (1 .-μ*a*abs.(k).^2).*(F₂.^2)
		end

		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end

		
		η = zeros(Complex{Float64}, mesh.N)
		h = zeros(Complex{Float64}, mesh.N)
		v = zeros(Complex{Float64}, mesh.N)
		fftη = zeros(Complex{Float64}, mesh.N)
		fftv = zeros(Complex{Float64}, mesh.N)


		# Evolution equations are ∂t U = f(U)
		function f!(U, currentTime)
			fftv .= F₂.*U[:,2]
			v .= real(ifft(fftv))

			fftη .= U[:,1]
		   	η .= real(ifft(U[:,1]))
			h .= 1 .+ ϵ .* η

			if noiseType == "Waves"
				stochVisc = sum(phiArray[:,2].^2) * smoothingField2
				fftus = ϵ .^2 ./ 2 .* ∂ₓ .* complexRounder.(fft(stochVisc))
				us = realRounder.(real.(ifft(fftus)))

				U[:,1] .= complexRounder.(- ∂ₓ.*(fftv.+ϵ*Π⅔.*fft(η.*v)) .+ ∂ₓ .* (fftus.+ϵ* Π⅔ .*fft(η.*us)))
				U[:,2] .= complexRounder.(- ∂ₓ.*(fftη.+ϵ/2*Π⅔.*fft(v.^2) ) .+ ϵ .* Π⅔ .* fft(us .* real.(ifft(∂ₓ .* fftv))))

				if !(isTrueStrato) && SGN
					@error "Type of simulation not supported: the Itô-Stratonovitch correction is not tractable in the " + label +" model (too complex).\
					Change parameter *isTrueStrato* to *true*, or remove it."
				elseif !(isTrueStrato)
					U[:,1] .+= ϵ .^ 2 ./ 2 .* ∂ₓ .* fft(stochVisc .* real.(ifft(∂ₓ .* fftη)))
					U[:,2] .+= ϵ .^ 2 ./ 2 .* ∂ₓ .* fft(stochVisc .* real.(ifft(∂ₓ .* fftv)))
				end
			else
				tablePhi = zeros(Complex{Float64}, (mesh.N, n_phi))
				stochVisc = zeros(Complex{Float64}, mesh.N)
				dG_IS = zeros(Complex{Float64}, mesh.N)

				for index in 1:mesh.N
					tablePhi[index,:] .= phi(x[index],currentTime)[:]
				end
				for k in 1:n_phi
					stochVisc .+= tablePhi[:,k] .^ 2
				end

				v .= real.(ifft(fftv))
				fftus = ϵ ./ 2 .* ∂ₓ .* fft(stochVisc)
				us = real.(ifft(fftus))

				U[:,1] .= -∂ₓ.*(fftv.+ϵ*Π⅔.*fft(η.*v)) .+ ∂ₓ .* (fftus.+ϵ*Π⅔ .*fft(η.*us))
				U[:,2] .= -∂ₓ.*(fftη.+ϵ/2*Π⅔.*fft(v .^2) ) .+ ϵ .* Π⅔ .* fft(us .* real.(ifft(∂ₓ .* fftv)))
				
			end
			
			U[abs.(U).< ktol ].=0
		end



		# Build raw data from physical data.
		# Discrete Fourier transform with, possibly, dealiasing and Krasny filter.
		function mapto(data::InitialData)
			U = [Π⅔ .* fft(data.η(x)) Π⅔ .*fft(data.v(x))]
			U[abs.(U).< ktol ].=0
			return U
		end

		# Reconstruct physical variables from raw data
		# Return `(η,v,x)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `x` is the vector of collocation points
		function mapfro(U)
			real(ifft(U[:,1])),real(ifft(U[:,2])),mesh.x
		end

		stoch_incr_U = zeros(Complex{Float64},(N, 2, 2*n_phi))

		function stoch_incr!(U,currentTime)
			
			fftη .= U[:,1]
			fftv .= F₂.*U[:,2]
			η .= real(ifft(U[:,1]))

			if noiseType == "Waves"

				for k in 1:n_phi
					temp_fftSigma = zeros(N)
					temp_fftSigma[1+round(Int, phiArray[k,1]/dk)] = phiArray[k,2]

					temp_sigma1 = realRounder.(N*(smoothingField .* real.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )
					temp_sigma2 = realRounder.(N*(smoothingField .* imag.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )

					stoch_incr_U[:,1,2*k-1] .+= - ∂ₓ .* (additiveNoise * fft(temp_sigma1) .+ ϵ * Π⅔ .* fft(η.*temp_sigma1))
					stoch_incr_U[:,2,2*k-1] .+= - ϵ .* Π⅔ .* fft(temp_sigma1 .* real.(ifft(∂ₓ .* fftv)))

					stoch_incr_U[:,1,2*k] .+= - ∂ₓ .* (additiveNoise * fft(temp_sigma2) .+ ϵ * Π⅔ .* fft(η.*temp_sigma2))
					stoch_incr_U[:,2,2*k] .+= - ϵ .* Π⅔ .* fft(temp_sigma2 .* real.(ifft(∂ₓ .* fftv)))
					
				end

			else
				tablePhi = zeros(Complex{Float64}, (mesh.N, n_phi))

				for i in 1:mesh.N
					tablePhi[i,:] .= phi(x[i],currentTime)[:]
				end

				for k in 1:n_phi
					stoch_incr_U[:,1,k] .+= - ϵ .* Π⅔ .* ∂ₓ .* (fft(tablePhi[:,k]).+ϵ*Π⅔.*fft(η.*tablePhi[:,k]))
					stoch_incr_U[:,2,k] .+= - ϵ .* Π⅔ .* fft( tablePhi[:,k] .* ∂ₓ .* ifft(fftv))
					if SGN
						dG = ϵ .* tablePhi[:,k] .* real.(ifft(∂ₓ .* ∂ₓ .* fftv)) .- ϵ .* real.(ifft(∂ₓ .* fft(tablePhi[:,k]))) .* real.(ifft( ∂ₓ .* fftv))
						stoch_incr_U[:,2,k] .+= ϵ .* μ .* Π⅔ .* fft( 1 ./ h .* real.(ifft( ∂ₓ .* fft(h.^3 ./3 .* dG))))
					end
				end
			end
			stoch_incr_U[abs.(stoch_incr_U).< ktol ].=0
			return stoch_incr_U
		end

		new(label, f!, mapto, mapfro, info, stoch_incr!, n_phi, isTrueStrato, smoothingField)
	end
end
