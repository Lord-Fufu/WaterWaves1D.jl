
using FFTW

function maxEta( pbm :: Problem)
	eta = real(ifft(last(pbm.data.U)))[:,1]
    l= length(eta)
    M = 0
    for k in 1:l
        M = max(eta[k], M)
    end
    return M
end

#=list_index=[]
for k in 1:n_sample
    if maxEta( pbmList_glob[k]) < 1
        push!(list_index,k)
    end
end

pbmList = pbmList_glob[list_index]=#

pbmList = pbmList_glob

#Compute the mean and the std of the problems (see utils-stat.jl in WaterWaves1D)
avg_pbm, std_pbm = meanAndStd(pbmList)
print("Mean and STD done.\n")

#Compute the maximal and the minimal problems (see utils-stat.jl in WaterWaves1D)
#max_pbm, min_pbm = maxAndMin(pbmList)
#print("Max and min done.\n")




#Build problems mean +- n*std
avgPlus2Std = avg_pbm + 2.0*std_pbm
avgMinus2Std = avg_pbm - 2.0*std_pbm
avgPlus3Std = avg_pbm + 3.0*std_pbm
avgMinus3Std = avg_pbm - 3.0*std_pbm

avgPlus2Std.label = "Average + 2 STD"
avgMinus2Std.label = "Average - 2 STD"
avgPlus3Std.label = "Average + 3 STD"
avgMinus3Std.label = "Average - 3 STD"

avg_pbm.label = "Average"
std_pbm.label = "STD"

print("Data to plot done.\n")
