#' Load model from cache
#'
#' Check if a model is cached and load it if it is, compute it otherwise
#' @export
checkModel <- function(fname, model_links, basal_file, data.stimulation, data.variation="", nb_cores=1, inits=1000, perform_plots=F, precorrelate=T, method="geneticlhs", unused_perturbations=c(), unused_readouts=c(), MIN_CV=0.1, DEFAULT_CV=0.3, model_name="default", data_space="linear") {

    # Check if a midas file exists
    if (fname %in% list.files(knitr::opts_chunk$get()$cache.path)) {
        model = suppressMessages( rebuildModel(paste0(knitr::opts_chunk$get()$cache.path, fname), data.stimulation, data.variation) )
    } else {
        model = createModel(model_links, basal_file, data.stimulation, data.variation, nb_cores, inits, perform_plots, precorrelate, method, unused_perturbations, unused_readouts, MIN_CV, DEFAULT_CV, model_name, data_space=data_space)
        exportModel(model, paste0(knitr::opts_chunk$get()$cache.path, fname))
    }
    print(paste("Best fit:", signif(model$bestfit, 5), ", Score=", signif(model$bestfitscore, 2)))
    return(model)
}

#' Check if a ModelSet exists or create it
#' @export
checkModelSet <- function(fname, model_links, basal_file, data.stimulation, data.variation=c(), nb_cores=1, inits=1000, perform_plots=F, method="geneticlhs", unused_perturbations=c(), unused_readouts=c(), MIN_CV=0.1, DEFAULT_CV=0.3, model_name="default", data_space="linear") {
    load_files = paste0(fname, "_", names(data.stimulation), ".mra")
    # Check if a midas file exists
    if (all(load_files %in% list.files(knitr::opts_chunk$get()$cache.path))) {
        model = suppressMessages( rebuildModelSet(paste0(knitr::opts_chunk$get()$cache.path, load_files), data.stimulation, data.variation) )
    } else {
        model = createModelSet(model_links, basal_file, data.stimulation, data.variation, nb_cores, inits, perform_plots, method, unused_perturbations, unused_readouts, MIN_CV, DEFAULT_CV, model_name, data_space=data_space)
        subm = extractSubmodels(model)
        for (ii in 1:length(subm$models)) {
            exportModel(subm$models[[ii]], paste0(knitr::opts_chunk$get()$cache.path, load_files[ii]))
        }
    }
    print(paste("Best fit:", signif(model$bestfit, 5), ", Score=", signif(model$bestfitscore, 2)))
    return(model)
}

#' Perform variable parameters analysis or load a variable model if it exists
#' @export
checkVariableParameters <- function(vmodel_name, original_modelset, data.stimulation, data.variation=c(), nb_cores=0, max_iterations=0, nb_samples=100, accuracy=0.95, method="geneticlhs", notVariable=c()) {
    load_files = paste0(vmodel_name, "_", names(data.stimulation), ".mra")
    # Check if a midas file exists
    if (all(load_files %in% list.files(knitr::opts_chunk$get()$cache.path))) {
        model = suppressMessages( rebuildModelSet(paste0(knitr::opts_chunk$get()$cache.path, load_files), data.stimulation, data.variation) )
    } else {
        model = addVariableParameters(original_modelset, nb_cores, max_iterations, nb_samples, accuracy, method, notVariable)
        subm = extractSubmodels(model)
        for (ii in 1:length(subm$models)) {
            exportModel(subm$models[[ii]], paste0(knitr::opts_chunk$get()$cache.path, load_files[ii]))
        }
    }
    print(paste("Best fit:", signif(model$bestfit, 5), ", Score=", signif(model$bestfitscore, 2)))
    return(model)
}