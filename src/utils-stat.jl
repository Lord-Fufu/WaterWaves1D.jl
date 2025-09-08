


export meanAndStd,maxAndMin



function meanAndStd(pbmList :: Vector{Any})
	N = size(pbmList)[1]
	avg = (1/N) * pbmList[1]
	for k in 2:N
		avg = avg + (1/N) * pbmList[k]
	end
	
	var = (1/(N-1)) * (pbmList[1] - avg)^2
	for k in 2:N
    	var = var + (1/(N-1)) * (pbmList[k] - avg)^2
	end

	std = √(var)

	return (avg,std)
end

function maxAndMin(pbmList :: Vector{Any})
	N = size(pbmList)[1]

	maxP = pbmList[1]
	minP = pbmList[1]
	for k in 2:N
		maxP = max(maxP, pbmList[k])
		minP = min(minP, pbmList[k])
	end

	return (maxP,minP)
end
