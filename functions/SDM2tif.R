# run MaxEnt and output tif
  # predictors: environment factor 
    # (predict and project layer, predict layer must be the first)
  # data: species records
  # path.sdm: output directory
  # maxent_args: maxent arguments

SDM2tif <- 
  function(predictors, 
           data,
           bias = NULL,
           path.sdm, 
           maxent_args){
    #-- run MaxEt
    SDM <- 
      maxent(predictors[[1]], data,
             path = path.sdm,
             a = bias,
             args = maxent_args,
             silent = TRUE)
    
    #-- export as tif
    lapply(1:length(predictors),
           function(x){
             ## 2010.stack
             SDM.stack <- 
               lapply(1:length(SDM@models),
                      function(i)
                        predict(predictors[[x]], SDM@models[[i]])) %>% 
               stack()
             writeRaster(SDM.stack, 
                         sprintf("%s/model_predict_stack_%s.tif",
                                 path.sdm, x))
             ## 2010.avg
             SDM.avg <- calc(SDM.stack, fun = mean)
             writeRaster(SDM.avg, 
                         sprintf("%s/model_predict_avg_%s.tif",
                                 path.sdm, x))
           })
    return(SDM)
  }
