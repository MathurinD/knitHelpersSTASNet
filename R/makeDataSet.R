#' Make a dataset to generate a ModelSet with the same perturbations
#'
#' TODO fix to actually only have common ones
#' TODO add the possibility to add the non common as NA
#' @export
makeDataSet <- function(...) {
    data_sets = list(...)
    shared_conditions = Reduce(intersect, sapply(data_sets, colnames))
    condition_codes = sapply(data_sets, function(dts){ dts %>% select(matches("^TR")) %>% unite("code", matches("TR")) })
    data_sets = lapply(data_sets, function(dd){
	    		if ( dd %>% select(-shared_conditions) %>% ncol > 0 ) {
	    			return(dd %>% filter_at(vars(-shared_conditions), all_vars(.==0)) %>% select(shared_conditions) %>% as.data.frame)
			} else {
				return( dd )
			}
		})
    return(data_sets)
}
