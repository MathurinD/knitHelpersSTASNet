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

#' Load refited model from cache
#'
#' Check if a model exists or use the refitModel function to create it. Provides some helpful functions to change just a few parameters
#' @param fixed_parameters A list of parameters to fix in the form 'list("param_name"=param_value)'
#' @export
checkRefitModel <- function(fname, initial_model, fixed_parameters, nb_cores=1, inits=1000) {
    new_params = initial_model$parameters
    fixed_ids = c()
    for (nn in names(fixed_parameters)) {
        nn_id = which(getParametersNames(initial_model) == nn)
        if (length(nn_id) == 0) { stop(paste("Can't fix parameter. Unknown parameter name: ", nn)) }
        new_params[nn_id] = fixed_parameters[[nn]]
        fixed_ids = c(fixed_ids, nn_id)
    }

    # Check if a midas file exists
    if (fname %in% list.files(knitr::opts_chunk$get()$cache.path)) {
        model = suppressMessages( rebuildModel(paste0(knitr::opts_chunk$get()$cache.path, fname)) )
    } else {
        model = refitModel(initial_model, new_params, (1:length(new_params))[-fixed_ids], inits=inits, nb_cores=nb_cores)
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

#' Perform and cache model extension or load it from cache
#'
#' All parameters except 'fname' correspond to the ones in suggestExtension
#' @param fname Name for the saved file (without extension)
#' @param print Provided for convenient copy-paste from suggestExtension lines, has no control here as the results will always be cached in a file.
#' @param parallel Kept for backward compatibility, won't be used
#' @export
checkExtension <- function(fname, original_model, mc = 0, sample_range=c(10^(2:-1),0,-10^(-1:2)), padjust_method="bonferroni", parallel=TRUE, print = TRUE){
    fname = paste0(knitr::opts_chunk$get()$cache.path, "/", fname, ".tsv")
    if (file.exists(fname)) {
        return(read.table(fname, header=TRUE))
    } else {
        return(suggestExtension(original_model=original_model, mc=mc, sample_range=sample_range, padjust_method=padjust_method, print=TRUE, fname=fname))
    }
}
#' Manage profile likelihood caching
#'
#' Check cache for an existing save of the profile likelihood labelled 'fname'. Load it if it exists or perform profileLikelihood and save it if it doesn't.
#' All arguments except 'fname' correspond to the ones of profileLikelihood
#' @param fname Name for the saved file (without extension)
#' @export
checkProfileLikelihood <- function(fname, model_description, nb_points=10000, nb_cores=1) {
    fname = paste0(knitr::opts_chunk$get()$cache.path, "/", fname, ".rds")
    if (file.exists(fname)) {
        return(importProfiles(fname))
    } else {
        pl = profileLikelihood(model_description, nb_points, nb_cores)
        exportProfiles(pl, fname)
        return(pl)
    }
}
