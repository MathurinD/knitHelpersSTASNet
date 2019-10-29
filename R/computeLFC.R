#' Compute log-fold changes from a MIDAS file
#' @export
computeLFC <- function(full_datas) {
        control_selection = grep("^c$|control", rownames(full_datas))
        control=full_datas[control_selection,,drop=FALSE]
        if (!is.null(nrow(control)) && nrow(control) > 0) {
            control_line = colMeans(control, na.rm=TRUE)
        } else {
            control_line = rep(1, ncol(full_datas))
        }
        blank_selection = which(rownames(full_datas)=="blank")
        blank=full_datas[blank_selection,]
        if (length(blank_selection) > 0 || length(control_selection) > 0) {
            datas=full_datas[-c( blank_selection, control_selection ),]
        } else {
            datas = full_datas
        }
        control = matrix(rep(control_line, nrow(datas)), nrow=nrow(datas), byrow=T)
        rep_variation = aggregate(full_datas, by=list(rownames(full_datas)), sd)
        rownames(rep_variation) = rep_variation[,1]
        rep_variation = rep_variation[,-1]
        rep_mean = aggregate(full_datas, by=list(rownames(full_datas)), mean)[,-1]
        return(log2(datas/control))
}

#' Plots to compare replicates
#'
#' Compare replicates with readouts as columns and treatments as rownames
#' @param rep1 Data for replicate 1
#' @param rep2 Data for replicate 2
#' @param rep1_name Name of replicate 1
#' @param rep2_name Name of replicate 2
#' @export
compare_replicates <- function(rep1, rep2, rep1_name="", rep2_name="") {
    rownames(rep1) = make.names(rownames(rep1))
    rownames(rep2) = make.names(rownames(rep2))
    shared_names = rownames(rep1)[rownames(rep1)%in%rownames(rep2)]
    for (rr in colnames(rep1)) {
        if (rr %in% colnames(rep2)) {
            plot(rep1[shared_names,rr], rep2[shared_names,rr], xlab=rep1_name, ylab=rep2_name, main=paste0("Replicates correlation for ", rr))
            #plot(log2(rep1[shared_names,rr]), log2(rep2[shared_names,rr]), xlab=rep1_name, ylab=rep2_name, main=paste0("Replicates correlation for ", rr))
            lines(-15:15, -15:15)
        }
    }
}
