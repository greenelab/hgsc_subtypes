setClass(Class = "DoppelGang",
         representation(fullresults = "list", summaryresults="data.frame", inputargs="list"))

setGeneric("print")
setGeneric("summary")
setGeneric("plot")
#setGeneric("show")

setMethod("summary", signature(object="DoppelGang"),
          function(object) object@summaryresults)

setMethod("show", signature(object="DoppelGang"),
          function(object){
              cat ("S4 object of class:" , class ( object ) , "\n")
              cat ("Number of potential doppelgangers:" , nrow(object@summaryresults), ": ",
                   sum(object@summaryresults$expr.doppel), "expression, ",
                   sum(object@summaryresults$pheno.doppel), "phenotype, ",
                   sum(object@summaryresults$smokinggun.doppel), "smoking gun. \n")
              cat (" \n ")
              cat ("Use summary(object) to obtain a data.frame of potential doppelgangrs.", "\n")
              cat (" \n ")
          })

setMethod("print", signature(x="DoppelGang"),
          function(x) print(x@summaryresults))


setMethod("plot", signature(x="DoppelGang"),
          function(x, skip.no.doppels=FALSE, plot.pair=NULL, ...){
              if(!is.null(plot.pair)){
                  if(is.character(plot.pair) & length(plot.pair) == 2){
                      study.names <- unique(unlist(strsplit(names(x@fullresults), split=x@inputargs$separator)))
                      if(!all(plot.pair %in% study.names))
                          stop("One or both of plot.pair do not match names(esets)")
                      iplot <- which(names(x@fullresults) %in% c(paste(plot.pair, collapse=":"), paste(rev(plot.pair), collapse=":")))
                  }else{
                      stop("plot.pair must be a character vector of length two, containing the names of two datasets seen in names(esets)")
                  }
              }else{
                  iplot <- 1:length(x@fullresults)
              }
              op <- par()
              for (i in iplot){
                  cors <- x@fullresults[[i]]$correlations
                  cors <- na.omit(as.numeric(cors))
                  expr.doppels <- x@fullresults[[i]]$expr.doppels$outlierFinder.res
                  expr.doppels <- expr.doppels[expr.doppels$doppel, ]
                  if(skip.no.doppels & (nrow(expr.doppels) == 0 | is.null(expr.doppels)))
                      next
                  stfit <- x@fullresults[[i]]$expr.doppels$stfit
                  hist(cors, main = paste(names(x@fullresults)[i], stfit$algorithm$message, sep="\n"),
                       freq=TRUE, xlab = "Pairwise Correlations", breaks="FD", ...)
                  abline(v=expr.doppels$similarity, col="red", lw=0.5)
                  xvals <- seq(min(-0.999, par("usr")[1]), max(0.999, par("usr")[2]), length.out=100)
                  transFun <- x@inputargs$outlierFinder.expr.args$transFun
                  suppressWarnings(xvals <- xvals[!is.na(transFun(xvals))])
                  yvals <- diff(pst(transFun(xvals), location=stfit$dp["location"], scale=stfit$dp["scale"], shape=stfit$dp["shape"], df=stfit$dp["df"])) * length(cors)
                  lines(xvals[-1], yvals)
                  par(ask = prod(par("mfcol")) < length(x@fullresults) && dev.interactive())
              }
              suppressWarnings(par(op))
          })
