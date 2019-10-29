#' Plot the values in a plate from a MIDAS file with ID:well column
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
