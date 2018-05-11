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
#par[1] = wavelength
#par[2] = phase
#par[3] = scale
#par[4] = offset
cosineFunction = function(par, x){
    x = x-par[2]
    return(
        (abs(par[1]*x) < pi) * (par[4]+par[3]*(1+cos(par[1]*x))) +
        (abs(par[1]*x) >= pi) * par[4]
    )
}

cosineTuningCurveMSE = function(par, data){
    estimate = cosineFunction(par,data$angle_bins)
    return(sum((estimate-data$firing_rate)^2))
}

addGaussianEstimate = function(cells){

    cells$gaussian_estimate = 0
    cells$gaussian_RMSE = 0
    for(cellname in unique(cells$cellname)){
        
        cell = cells[cells$cellname == cellname,]
        initial_values = c(mean=cell$angle_bins[which.max(cell$firing_rate)], sd=50, scale=10, offset=0)
        
        result = optim(initial_values, gaussianTuningCurveMSE, data=cell, method="BFGS", control=list(maxit=1000))
        if(result$convergence != 0){
            print("Did not converge")
        }
        
        cells[cells$cellname == cellname,]$gaussian_estimate = gaussianFunction(result$par, cell$angle_bins)
        cells[cells$cellname == cellname,]$gaussian_RMSE = sqrt(gaussianTuningCurveMSE(result$par, cell))
    }
    return(cells)
}

addCosineEstimate = function(cells){
    cells$cosine_estimate = 0
    cells$cosine_RMSE = 0
    for(cellname in unique(cells$cellname)){
        
        cell = cells[cells$cellname == cellname,]
        initial_values = c(wavelength=1/20, phase=cell$angle_bins[which.max(cell$firing_rate)], scale=max(cell$firing_rate)/2, offset=0)
        
        result = optim(initial_values, cosineTuningCurveMSE, data=cell, method="BFGS", control=list(maxit=1000))
        if(result$convergence != 0){
            print("Did not converge")
        }
        cells[cells$cellname == cellname,]$cosine_estimate = cosineFunction(result$par, cell$angle_bins)
        cells[cells$cellname == cellname,]$cosine_RMSE = sqrt(cosineTuningCurveMSE(result$par, cell))
    }
    return(cells)
}
