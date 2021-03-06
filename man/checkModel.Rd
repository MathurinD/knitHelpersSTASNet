% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caching.R
\name{checkModel}
\alias{checkModel}
\title{Load model from cache}
\usage{
checkModel(
  fname,
  model_links,
  basal_file,
  data.stimulation,
  data.variation = "",
  nb_cores = 1,
  inits = 1000,
  perform_plots = F,
  precorrelate = T,
  method = "geneticlhs",
  unused_perturbations = c(),
  unused_readouts = c(),
  MIN_CV = 0.1,
  DEFAULT_CV = 0.3,
  model_name = "default",
  data_space = "linear"
)
}
\description{
Check if a model is cached and load it if it is, compute it otherwise
}
