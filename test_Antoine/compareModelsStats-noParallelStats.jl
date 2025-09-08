

using CSV, DelimitedFiles
using FFTW, Statistics

n_begin = 1
n_end = 130

matrixOfTrajectories = zeros(Float64, 2048, n_end - n_begin + 1)
for p in 1:n_end
    #data_from_file = CSV.read(".//Documents//__plot Julia//SGN-files//traj" * string(p) * ".csv", DataFrame,
                                                                        #header=false, types=Complex{Float64})
    data_from_file = readdlm(".//Documents//__plot Julia//trajectoriesForStats//__other-SGN-files-p3-mediumNoise//traj" * string(p) * ".csv")
    
    arr_spectral=zeros(ComplexF64, 2048)
    arr_spectral_imag=zeros(ComplexF64, 2048)
    arr_sign=zeros(ComplexF64, 2048)

    arr_sign .= 2 .* (data_from_file[:,2].=="+") .- 1
    for p in 1:2048
        len = length(data_from_file[p,3])
        arr_spectral_imag[p] = parse(Float64, String(data_from_file[p,3][1:(len-2)]))
    end

    arr_spectral.= data_from_file[:,1] .+ arr_sign .* arr_spectral_imag .* 1im
    arr_for_plot = real.(ifft(arr_spectral))
    matrixOfTrajectories[:,p] .= arr_for_plot 
end

meanTraj = mean(matrixOfTrajectories, dims=2)
stdTraj  = std( matrixOfTrajectories, dims=2)




using Plots

T_plot = 5
index_to_plot = 11

_, _, meshToPlot = det_problem.model.mapfro(det_problem.data.U[index_to_plot])

plot(meshToPlot, meanTraj, ribbon = 3 .* stdTraj,
            fillalpha = 0.5, c=:orange, label = "3 STD Spread", legend = :topright,
            xlabel = "Horizontal domain (x)", ylabel = "Surface elevation (η)",
            titlefontsize=14,
            guidefontsize=14,
            tickfontsize=10,
            legendfontsize=10,
            xlim = (-20,20),
            ylim = (-0.4,1.0));
plot!(meshToPlot, meanTraj .+ 3 .* stdTraj, c=:orange, label = nothing);
plot!(meshToPlot, meanTraj .- 3 .* stdTraj, c=:orange, label = nothing);

plot!(meshToPlot, meanTraj, xlims=(-20,20), c=:red, label = "Average LU-SGN");
plot!(WW_problem, xlims=(-20,20), T=T_plot, c=:green, alpha = 0.5, label = "Water Waves")
plot!(det_problem, xlims=(-20,20), T=T_plot, c=:blue, alpha = 0.5, label = "SGN",
        xlabel = "Horizontal domain (x)", ylabel = "Surface elevation (η)")

savefig(".//Documents//__plot Julia//plotsForPaper//LUSGN-statistics.pdf")
