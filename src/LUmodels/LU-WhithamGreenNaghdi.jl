export LUWhithamGreenNaghdi

"""
    WhithamGreenNaghdi(param;kwargs)

Define an object of type `AbstractModel` in view of solving the initial-value problem for
the fully dispersive Green-Naghdi model proposed by [Duchêne, Israwi and Talhouk](https://doi.org/10.1137/130947064).

# Argument
`param` is of type `NamedTuple` and must contain
- dimensionless parameters `ϵ` (nonlinearity) and `μ` (dispersion);
- numerical parameters to construct the mesh of collocation points, if `mesh` is not provided as a keyword argument.

## Optional keyword arguments
- `SGN`: if `true` (default is `false`), compute the Serre-Green-Naghdi (SGN) instead of Whitham-Green-Naghdi (WGN) system (see `SerreGreenNaghdi(param;kwargs)`);
- `mesh`: the mesh of collocation points. By default, `mesh = Mesh(param)`;
- `iterative`: solve the elliptic problem through GMRES if `true`, LU decomposition if `false` (default is `true`);
- `precond`: use a (left) preconditioner for GMRES if `true` (default), choose `precond` as the preconditioner if provided;
- `gtol`: relative tolerance of the GMRES algorithm (default is `1e-14`);
- `restart`: the corresponding option of the GMRES algorithm (default is `100`);
- `maxiter`: the corresponding option of GMRES (default is `nothing`);
- `ktol`: tolerance of the Krasny filter (default is `0`, i.e. no filtering);
- `dealias`: dealiasing with Orlicz rule `1-dealias/(dealias+2)` (default is `0`, i.e. no dealiasing);
- `label`: a label for future references (default is `"Whitham-Green-Naghdi"`);

# Return values
Generate necessary ingredients for solving an initial-value problem via `solve!`:
1. a function `WhithamGreenNaghdi.f!` to be called in explicit time-integration solvers;
2. a function `WhithamGreenNaghdi.mapto` which from `(η,v)` of type `InitialData` provides the raw data matrix on which computations are to be executed;
3. a function `WhithamGreenNaghdi.mapfro` which from such data matrix returns the Tuple of real vectors `(η,v,x)`, where
	- `η` is the values of surface deformation at collocation points `x`;
	- `v` is the derivative of the trace of the velocity potential at `x`;
4. additionally, a handy function `WhithamGreenNaghdi.mapfrofull` which from data matrix returns the Tuple of real vectors `(η,v,u)`, where
    - `u` corresponds to the layer-averaged velocity.

"""
mutable struct LUWhithamGreenNaghdi <: AbstractModel

	label   :: String
	f!		:: Function
	mapto	:: Function
	mapfro	:: Function
	mapfrofull	:: Function
	info    :: String
	stoch_incr! :: Function
	n_phi   :: Int
	isTrueStrato :: Bool
	noiseType :: String


    function LUWhithamGreenNaghdi(param::NamedTuple; SGN = false,
								dealias	= 0,
								ktol	= 0,
								iterate	= true,
								gtol	= 1e-14,
								precond	= true,
								restart	= nothing,
								maxiter	= nothing,
								label	= nothing
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

		if isnothing(maxiter) maxiter = mesh.N end
		if isnothing(restart) restart = min(20,mesh.N) end
		if isnothing(label)
			if SGN == true
				label = "LU-Serre-Green-Naghdi"
			else
				label = "LU-Whitham-Green-Naghdi"
			end
		end

		if !isTrueStrato
			@error "Type of simulation not supported: the Itô-Stratonovitch correction is not tractable in the " + label +" model (too complex).\
			Change parameter *isTrueStrato* to *true*, or remove it."
		end

		if in(:phiArray,keys(param))
			phiArray = param.phiArray
			smoother = param.smoother
			smoothingField = smoother.(mesh.x)
			smoothingField2 = smoothingField.^2
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
		info = "$label model.\n"
		info *= "├─Shallowness parameter μ=$μ, nonlinearity parameter ϵ=$ϵ.\n"
		if dealias == 0
			info *= "├─No dealiasing. "
		else
			info *= "├─Dealiasing with Orszag's rule adapted to power $(dealias + 1) nonlinearity. "
		end
		if ktol == 0
			info *= "No Krasny filter. "
		else
			info *= "Krasny filter with tolerance $ktol."
		end
		if iterate == true
			if precond == false out="out" else out="" end
			info *= "\n└─Elliptic problem solved with GMRES method with$out preconditioning, \
			tolerance $gtol, maximal number of iterations $maxiter, restart after $restart iterations \
			(consider `iterate=false` for non-iterative method). "
		else
			info *= "\n└─Elliptic problem solved with standard LU factorization \
			(consider `iterate=true` for faster results). "
		end
		info *= "\nDiscretized with $(mesh.N) collocation points on [$(mesh.xmin), $(mesh.xmax)]."

		# Pre-allocate useful data
		k = mesh.k
		x = mesh.x
		x₀ = mesh.x[1]
		dk = mesh.dk
		N  = mesh.N

		∂ₓ	=  1im * k
        if SGN == true
			F₁ = 1 ./(1 .+ μ/3*k.^2)
			F₀ = sqrt(μ)*∂ₓ
	    else
			F₁ 	= tanh.(sqrt(μ)*abs.(k))./(sqrt(μ)*abs.(k))
			F₁[1] 	= 1
			F₀ = 1im * sqrt.(3*(1 ./F₁ .- 1)).*sign.(k)
		end
		if precond == true
			#Precond = Diagonal( 1 .+ μ/3*k.^2 )
			Precond = Diagonal( 1 ./  F₁ )
		elseif precond == false
			Precond = Diagonal( ones(size(k)) )
		else
			Precond = precond
		end
		if dealias == 0
			Π⅔ 	= ones(size(k)) # no dealiasing (Π⅔=Id)
		else
			K = (mesh.kmax-mesh.kmin)/(2+dealias)
			Π⅔ 	= abs.(k) .<= K # Dealiasing low-pass filter
		end
		FFT = exp.(-1im*k*(x.-x₀)');
        IFFT = exp.(1im*k*(x.-x₀)')/length(x);
		M₀ = IFFT * Diagonal( F₀ .* Π⅔)* FFT
		IFFTF₀ = IFFT * Diagonal( F₀ .* Π⅔)
        Id = Diagonal(ones(size(x)));
		h = zeros(Complex{Float64}, mesh.N)
		u, η, fftv, fftη, fftu, hdu = (similar(h),).*ones(6)
		L = similar(FFT)
		dG_IS = zeros(Complex{Float64}, mesh.N)		


		# Evolution equations are ∂t U = f(U)
		function f!(U, currentTime)
			fftη .= U[:,1]
			h .= 1 .+ ϵ*ifft(fftη)
			fftv .= U[:,2]

			if iterate == false
				L .= Id - 1/3 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFTF₀
				fftu .= L \ fftv
			elseif iterate == true
		        function LL(hatu)
		            hatu- 1/3 *Π⅔.*fft( 1 ./h .* ifft( F₀ .* Π⅔.*fft( h.^3 .* ifft( F₀ .* Π⅔.* hatu ) ) ) )
				end
				fftu .= gmres( LinearMap(LL, length(h); issymmetric=false, ismutating=false) , fftv ;
						restart = restart, maxiter = maxiter, Pl = Precond, reltol = gtol )
			end


			if noiseType == "Waves"
				u .= ifft(fftu)
				hdu .= h .* ifft(Π⅔.*F₀.*fftu)

				stochVisc = sum(phiArray[:,2].^2) * smoothingField2
				fftus = 1 ./ 2 .* ∂ₓ .* fft(stochVisc)
				us = ifft(fftus)

				dG_IS  .= -us .* ifft(∂ₓ .* ∂ₓ .* fftv) .+ ifft(∂ₓ .* fftus) .* ifft(∂ₓ .* fftv)		
				U[:,1] .= -∂ₓ.*Π⅔.*(fftu .+ ϵ * fft(ifft(fftη) .* u)) .+ ϵ .* ∂ₓ .* fft(h .* us )
				U[:,2] .= -∂ₓ.*Π⅔.*(fftη .+ ϵ * fft( u.*ifft(fftv) .- 1/2 * u.^2 .- 1/2 * hdu.^2 ) )
						 .+ ϵ .* fft(us .* ifft(∂ₓ .* fftv)) .+ ϵ^2 * μ ./ h .* ∂ₓ .* fft(h .^3 ./3 .* dG_IS)

				U[abs.(U).< ktol ].=0
				
			else
				tablePhi = zeros(Complex{Float64}, (mesh.N, n_phi))
				stochVisc = zeros(Complex{Float64}, mesh.N)
				us = zeros(Complex{Float64}, mesh.N)
				fftus = zeros(Complex{Float64}, mesh.N)
				dG_IS = zeros(Complex{Float64}, mesh.N)

				for index in 1:mesh.N
					tablePhi[index,:] .= phi(x[index],currentTime)[:]
				end
				for k in 1:n_phi
					stochVisc .+= tablePhi[:,k] .^ 2
				end

				u .= ifft(fftu)
				hdu .= h .* ifft(Π⅔.*F₀.*fftu)

				fftus .= ϵ ./ 2 .* ∂ₓ .* fft(stochVisc)
				us .= ifft(fftus)

				dG_IS  .= -us .* ifft(∂ₓ .* ∂ₓ .* fftu) .+ ifft(∂ₓ .* fftus) .* ifft(∂ₓ .* fftu)		
				U[:,1] .= -∂ₓ.*Π⅔.*(fftu .+ ϵ * fft(ifft(fftη) .* u)) .+ ϵ .* ∂ₓ .* fft(h .* us )
				U[:,2] .= -∂ₓ.*Π⅔.*(fftη .+ ϵ * fft( u.*ifft(fftv) .- 1/2 * u.^2 .- 1/2 * hdu.^2 ) )
													.+ ϵ .* fft(us .* ifft(∂ₓ .* fftu)) .+ ϵ^2 * μ ./ h .* ∂ₓ .* fft(h .^3 ./3 .* dG_IS)
				U[abs.(U).< ktol ].=0
			end
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
		# Return `(η,v,u)`, where
		# - `η` is the surface deformation;
		# - `v` is the derivative of the trace of the velocity potential;
		# - `u` corresponds to the layer-averaged velocity.
		# Inverse Fourier transform and take the real part, plus solves the costly elliptic problem for `u`.
		function mapfrofull(U)
				fftη .= U[:,1]
			   	h .= 1 .+ ϵ*ifft(fftη)
				L .= Id - 1/3 * Diagonal(Π⅔) * FFT * Diagonal( 1 ./h ) * M₀ * Diagonal( h.^3 ) * IFFTF₀
				real(ifft(U[:,1])),real(ifft(U[:,2])),real(ifft(L \ U[:,2]))
		end
		

		#=function stoch_incr!(U,currentTime)

			if !isTrueStrato
				@error "Methods other than 'true Strato' (i.e Euler-Heun) are not supported for this model. Change parameter 'isTrueStrato' to 'true'."
			end

			fftη .= U[:,1]
			fftv .= U[:,2]
			η .= ifft(fftη)
			h .= 1 .+ ϵ*ifft(fftη)

			stoch_incr_U = zeros(Complex{Float64},(mesh.N, 2, 2*n_phi))


			if noiseType == "Waves"

				for k in 1:n_phi

					temp_fftSigma = zeros(mesh.N)
					temp_fftSigma[1+round(Int, phiArray[k,1]/dk)] = phiArray[k,2]

					temp_sigma1 = realRounder.(mesh.N*(real.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )
					temp_sigma2 = realRounder.(mesh.N*(smoothingField .* imag.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )

					stoch_incr_U[:,1,2*k-1] .+= - ϵ .* ∂ₓ .* (additiveNoise * fft(temp_sigma1) .+ ϵ * Π⅔ .* fft(η.*temp_sigma1))
					stoch_incr_U[:,2,2*k-1] .+= - ϵ .* Π⅔ .* fft(temp_sigma1 .* real.(ifft(∂ₓ .* fftv)))

					stoch_incr_U[:,1,2*k] .+= - ϵ .* ∂ₓ .* (additiveNoise * fft(temp_sigma2) .+ ϵ * Π⅔ .* fft(η.*temp_sigma2))
					stoch_incr_U[:,2,2*k] .+= - ϵ .* Π⅔ .* fft(temp_sigma2 .* real.(ifft(∂ₓ .* fftv)))
					
					dG1 = temp_sigma1 .* real.(ifft(∂ₓ .* ∂ₓ .* fftv)) .- real.(ifft(∂ₓ .* fft(temp_sigma1))) .* real.(ifft( ∂ₓ .* fftv))
					dG2 = temp_sigma2 .* real.(ifft(∂ₓ .* ∂ₓ .* fftv)) .-  real.(ifft(∂ₓ .* fft(temp_sigma2))) .* real.(ifft( ∂ₓ .* fftv))
					stoch_incr_U[:,2,2*k-1] .+= ϵ^2 * μ .* Π⅔ .* fft( 1 ./ h .* real.(ifft( ∂ₓ .* fft(h.^3 ./3 .* dG1))))
					stoch_incr_U[:,2,2*k]   .+= ϵ^2 * μ .* Π⅔ .* fft( 1 ./ h .* real.(ifft( ∂ₓ .* fft(h.^3 ./3 .* dG2))))
				end

			else
				tablePhi = zeros(Complex{Float64}, (mesh.N, n_phi))

				for i in 1:mesh.N
					tablePhi[i,:] .= phi(x[i],currentTime)[:]
				end

				for k in 1:n_phi
					stoch_incr_U[:,1,k] .+= - ϵ .* Π⅔ .*  fft( tablePhi[:,k] .* ifft(∂ₓ .* fftη))
					dG = tablePhi[:,k] .* ifft(∂ₓ .* ∂ₓ .* fftv) .- ifft(∂ₓ .* fft(tablePhi[:,k])) .* ifft( ∂ₓ .* fftv)
					stoch_incr_U[:,2,k] .+= - ϵ .* Π⅔ .* fft( tablePhi[:,k] .* ifft(∂ₓ .* fftv))
										   .+ ϵ^2 * μ .* Π⅔ .* fft( 1 ./ h .* ifft( ∂ₓ .* fft(h.^3 ./3 .* dG)))
				end
			end

			return stoch_incr_U
		end=#

		stoch_incr_U = zeros(Complex{Float64},(N, 2, 2*n_phi))

		function stoch_incr!(U,currentTime)
			
			fftη .= U[:,1]
			fftv .= U[:,2]
			η .= real(ifft(U[:,1]))

			if noiseType == "Waves"

				for k in 1:n_phi
					temp_fftSigma = zeros(N)
					temp_fftSigma[1+round(Int, phiArray[k,1]/dk)] = phiArray[k,2]

					temp_sigma1 = realRounder.(N*(smoothingField .* real.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )
					temp_sigma2 = realRounder.(N*(smoothingField .* imag.(ifft(temp_fftSigma)))) #* exp(-1im * phiArray[k,1] * currentTime ) )

					stoch_incr_U[:,1,2*k-1] .+= -  ∂ₓ .* (additiveNoise * fft(temp_sigma1) .+ ϵ * Π⅔ .* fft(η.*temp_sigma1))
					stoch_incr_U[:,2,2*k-1] .+= - ϵ .* Π⅔ .* fft(temp_sigma1 .* real.(ifft(∂ₓ .* fftv)))

					stoch_incr_U[:,1,2*k] .+= - ∂ₓ .* (additiveNoise * fft(temp_sigma2) .+ ϵ * Π⅔ .* fft(η.*temp_sigma2))
					stoch_incr_U[:,2,2*k] .+= - ϵ .* Π⅔ .* fft(temp_sigma2 .* real.(ifft(∂ₓ .* fftv)))

					dG1 = temp_sigma1 .* real.(ifft(∂ₓ .* ∂ₓ .* fftv)) .- real.(ifft(∂ₓ .* fft(temp_sigma1))) .* real.(ifft( ∂ₓ .* fftv))
					dG2 = temp_sigma2 .* real.(ifft(∂ₓ .* ∂ₓ .* fftv)) .-  real.(ifft(∂ₓ .* fft(temp_sigma2))) .* real.(ifft( ∂ₓ .* fftv))
					stoch_incr_U[:,2,2*k-1] .+= ϵ^2 * μ .* Π⅔ .* fft( 1 ./ h .* real.(ifft( ∂ₓ .* fft(h.^3 ./3 .* dG1))))
					stoch_incr_U[:,2,2*k]   .+= ϵ^2 * μ .* Π⅔ .* fft( 1 ./ h .* real.(ifft( ∂ₓ .* fft(h.^3 ./3 .* dG2))))					
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

        new(label, f!, mapto, mapfro, mapfrofull, info, stoch_incr!, n_phi, isTrueStrato )
    end
end
