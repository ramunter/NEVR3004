#Input
#par[1] = mean
#par[2] = sd
#par[3] = scale
#par[4] = offset
gaussianFunction = function(par, x){
    return(par[3]*exp(-((x-par[1])/par[2])^2) + par[4])
}

gaussianTuningCurveMSE = function(par, data){
    estimate = gaussianFunction(par,data$angle_bins)
    return(sum((estimate-data$firing_rate)^2))
}

#Input
#par[1] = amplitude
#par[2] = phase
#par[3] = scale
#par[4] = offset
cosineFunction = function(par, x){
    return(
        (abs(par[1]*(x-par[2])) < x) * par[4]+par[3]*(1+cos(par[1]*(x-par[2]))) +
        (abs(par[1]*(x-par[2])) >= x) * par[4]
    )
}

cosineTuningCurveMSE = function(par, data){
    estimate = cosineFunction(par,data$angle_bins)
    return(sum((estimate-data$firing_rate)^2))
}

addGaussianEstimate = function(cells){
    initial_values = c(mean=0, sd=100, scale=1, offset = 150)
    cells$gaussian_estimate = 0
    for(cellname in cells$cellname){
        cell = cells[cells$cellname == cellname,]
        result = optim(initial_values, gaussianTuningCurveMSE, data=cell, method="BFGS")
        cells[cells$cellname == cellname,]$gaussian_estimate = gaussianFunction(result$par, cell$angle_bins)
    }
    return(cells)
}

addCosineEstimate = function(cells){
    initial_values = c(amplitude=0, phase=100, scale=1, offset = 150)
    cells$cosine_estimate = 0
    for(cellname in cells$cellname){
        cell = cells[cells$cellname == cellname,]
        result = optim(initial_values, cosineTuningCurveMSE, data=cell, method="BFGS")
        cells[cells$cellname == cellname,]$cosine_estimate = cosineFunction(result$par, cell$angle_bins)
    }
    return(cells)
}
