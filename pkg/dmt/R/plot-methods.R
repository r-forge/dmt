plot.DependencyModel <- function(x, X, Y = NULL, 
		     ann.types = NULL, ann.cols = NULL, 
		     legend.x = 0, legend.y = 1, 
                     legend.xjust = 0, legend.yjust = 1, 
		     order=FALSE, 
		     cex.z = 0.6, cex.WX = 0.6, cex.WY = 0.6,...){

  model <- x
  
  # Check that both data sets are given for models from 2 data sets
  if (!is.null(getW(model)$X)){
    if (missing(X) || missing(Y)) {
      stop("Original data needed as 'X' and 'Y' arguments")
    }
  }
  # Check that data set is given for models from 1 data set
  else {
    if (missing(X)) {
      stop("Original data needed as 'X' argument")
    }
  }
  z <- z.effects(model,X,Y)[,1]
  W <- W.effects(model,X,Y)

  if (order){
    z.order <- order(z)
    z <- z[as.numeric(z.order)]
  }
  
  # Colors for different annotation types
  cols <- 'grey'
  if (!is.null(ann.types)) {
    if (length(ann.types) != ncol(X$data)) {
      warning("Length of ann.types doesn't match samples in data")
    }
    else {
      labels <- levels(ann.types) 
      types <- max(as.integer(ann.types),na.rm=TRUE)
      if (any(is.na(ann.types))) {
        types <- types + 1
	labels <- c(labels,"NA")
      }
      if (is.null(ann.cols)) {
        ann.cols <- gray(0:types / types)[1:types]
      }
      ann.int <- as.integer(ann.types)
      ann.int <- ifelse(is.na(ann.int),types,ann.int)
      cols <- ann.cols[ann.int]
      if (order) {
        cols <- cols[as.numeric(z.order)]
	ann.int <- ann.int[as.numeric(z.order)]
      }
    }
  }

  def.par <- par(no.readonly = TRUE) # save default

  # Outer margins and layout
  par(oma = c(0.5,0.5,2.5,0.5))

  # For models from 2 data sets
  if (!is.null(Y)){
    if (length(W$X) == 1){
      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(1,3),c(8,6))
    }
    else {
      layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE), c(2,2),c(8,6))
    }
  }
  # For models from 1 data set
  else {
    layout(matrix(c(1,2),2,1,byrow = TRUE), 1,c(8,6))
  }

  # Inner margins
  par(mar = c(3.1,4.1,3.1,0),cex.main = 1,cex.axis = 1)

  # Barplots for z
  barplot(z,main="Sample effects",col = cols,xlab = "Sample", ylab = "Weight",las = 2,cex.names=cex.z)  

  # Put legend 
  if(!is.null(ann.types)){
    legend(legend.x, legend.y, labels, cex=0.8, ann.cols, xjust=legend.xjust, yjust=legend.yjust);
  }

  par(mar = c(5.1,4.1,3.1,0.5),cex.main = 1,cex.axis = 1)
  # Barplots for Ws
  # For models from 2 data sets
  if (!is.null(Y)){
    if(length(W$X) == 1){
      barplot(t(W$X),main="Variable effects,\n(first data set)",ylab="Weight",las=2,cex.names=cex.WX)
    }
    else {
      barplot(t(W$X),main="Variable effects (first data set)",ylab="Weight",las=2,cex.names=cex.WX)
    }
    barplot(t(W$Y),main="Variable effects (second data set)",ylab="Weight",las=2,cex.names=cex.WY)
  }
  # For models from 1 data set
  else {
    barplot(t(W$total),main="Variable effects",ylab="Weight",las=2,cex.names=cex.WX)
  }

  #Title
  title <- paste(getModelMethod(model), "model")
  mtext(title, NORTH<-3, line=0, adj=0.5, cex=1.2, outer=TRUE)
  par(def.par)
}





plot.DependencyScreenModels <-
function (x, hilightGenes = NULL, showDensity = FALSE, 
	 showTop = 0, topName = FALSE,
         type = 'l', xlab = 'location', ylab = 'dependency score',
         main = paste('Dependency score'),
         pch = 20, cex = 0.75, tpch = 3, tcex = 1, ylim = NA, ...){

  models <- x
          
  #No printing without any models
  if(getModelNumbers(models) == 0) return()
        
  scores <- getScore(models)

  if(all(is.na(ylim))){
    ylim <- c(0,max(scores))
  }

  plot((locs/1e6), scores, type = type, xlab=xlab, 
       ylab=ylab, main=main, ylim = ylim,...)
        
  if(showTop > 0){
                
    d <- data.frame(scores,locs)
    d <- d[order(scores,decreasing=TRUE),]
                
    #Put vertical dashed line in the middle of showTop highest and showTop+1 highest scores
    limit = (d$score[showTop]+d$score[showTop+1])/2
    abline(h=limit,lty=2)
    topMods <- topModels(models,showTop) 
    #Draw points and add ranking number
    for(i in 1:showTop){
      points(d$locs[i]/1e6,d$scores[i],pch=tpch,cex = tcex) 
      if(topName)
        text(d$locs[i]/1e6,d$scores[i],getFeatureName(topMods[[i]]),pos=4,cex=tcex)
      else
        text(d$locs[i]/1e6,d$scores[i],as.character(i),pos=2)
    }
  }
        
  if (!is.null(hilightGenes)){

  #indices of cancer genes
  matches <- match(hilightGenes, featureNames)

  #print(matches)
  points(locs[matches]/1e6,scores[matches],pch=pch,cex=cex)

  }

  #lines to bottom to show gene density
  if(showDensity){
    heightCoefficient <- 100
    for(i in 1:length(locs)){
      lines( c((locs[i]/1e6) , (locs[i]/1e6)) , c(0,ylim[2]/heightCoefficient))
    }
  }
}
