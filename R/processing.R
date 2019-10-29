#' Compute log-fold change
#'@export
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

#' Plot values in plate layout
#' @export
plotPlateMIDAS <- function(midas_file) {
    midas_file = extractMIDAS(midas_file)
    wellX = substr(as.character(midas_file$ID.well), 1, 1)
    wellY = substr(as.character(midas_file$ID.well), 2, 10)
    readouts = colnames(midas_file)[grep("^DV", colnames(midas_file))]
    for (cc in readouts) {
        .GlobalEnv$plate = matrix( NA, nrow=8, ncol=12, dimnames=list(c("A", "B", "C", "D", "E", "F", "G", "H"), 1:12) ) # Empty plate
        dd = data.frame(wellX=wellX, wellY=wellY, Value=midas_file[,cc])
        rr=apply(dd, 1, function(xx){.GlobalEnv$plate[xx["wellX"], xx["wellY"]]=as.numeric(xx[3]) })
        plotHeatmap( log2(plate), main=paste0("Value across wells of ", gsub("DV.", "", cc)), lim=13, fixedRange=TRUE )
    }
}

