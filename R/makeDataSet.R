#' Make a dataset to generate a ModelSet with the same perturbations
#'
#' TODO fix to actually only have common ones
#' TODO add the possibility to add the non common as NA
#' @export
makeDataSet <- function(...) {
    data_sets = list(...)
    shared_conditions = Reduce(intersect, sapply(data_sets, colnames))
    data_sets = lapply(data_sets, function(dd){
	    		if ( dd %>% select(-shared_conditions) %>% ncol > 0 ) {
				to_keep = dd %>% select(-shared_conditions) %>% select(matches("TR")) %>% mutate(ID=1:nrow(.)) %>% filter_at(vars(-ID), all_vars(.==0)) %>% pull(ID)
				return( dd  %>% mutate(ID=1:nrow(.)) %>% filter(ID %in% to_keep) )
			} else {
				return( dd )
			}
		})
    return(data_sets)
}
