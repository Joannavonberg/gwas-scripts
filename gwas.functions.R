require(readr)

read.col <- function(dir = "./", fn, col.names, col.types){
  col.list <- list(col.types)
  names(col.list) <- col.names
  input <- paste0(dir, fn)
  return(read_table2(file = input, col_names = T, col_types = do.call(cols_only, col.list))[,names(col.list)])
}

stat.to.z <- function(stat, stat.type){
  z <- switch(stat.type,
              "Z" = stat,
              "P" = qnorm(stat/2),
              "CHISQ" = sqrt(stat)
  )
  return(z)
}

lam <- function(z, expected.stat){
  lambda = signif(median(z,na.rm=T)^2/expected.stat,4)
  return(lambda)
}

PlotQQHelper <- function(p, color, cex){
  p <- sort(p)
  lp <- length(p)
  lobs <- -log10(p)
  lexp <- -(log10(1:lp / (lp+1)))
  
  # plots all points with p < 0.05
  p.sig <- p < 0.05
  points(lexp[p.sig], lobs[p.sig], pch = 23, cex = cex, col = color, bg = color)
  
  # samples 5,000 points from p > 0.05 (to keep file size down)
  n=min(5000,sum(!p.sig))
  #ind <- sample((sum(p.sig)+1):lp, size = n, replace = F)
  ind <- seq(from = sum(p.sig)+1, to = lp, length.out = n)
  points(lexp[ind], lobs[ind], pch = 23, cex = cex, col = color, bg = color)
}

PlotQQ <- function(dir, infn, stratifyColName = NA, pvalColName = NA, 
                  stratifyColType = "n", pvalColType = "n", freqtomaf = T, 
                  colors = c("lightpink2", "purple", "deepskyblue1", "slateblue3", "olivedrab2"), 
                  fs = c(0, 0.001, 0.01, 0.05, 0.20, 0.50), obscol = rgb(255,79,0,maxColorValue=255), 
                  save = F, outfn = NA, expected.stat = 0.4549364){
  
  pval <- unlist(read.col(dir = dir, fn = infn, col.names = pvalColName, col.types = pvalColType))
  z <- stat.to.z(pval, "P")
  lambda <- lam(z, expected.stat)
  
  if(save){png(file = outfn, width = 8, height = 8, unit = "cm", pointsize = 4, res = 300)}
  
  plot(x = c(0,8), y = c(0,8), col = "gray25", lwd=3, type="l", xlab="Expected Distribution (-log10 of P value)", ylab="Observed Distribution (-log10 of P value)", xlim=c(0,8), ylim=c(0,8), las=1, xaxs="i", yaxs="i", bty="l", main=c(substitute(paste("QQ plot: ",lambda," = ", lam),list(lam = lambda)),expression()))
  PlotQQHelper(p = pval, color = obscol, cex = 0.4);
  if(!is.na(stratifyColName)){
    strat <- unlist(read.col(dir = dir, fn = infn, col.names = stratifyColName, col.types = stratifyColType))
    if(freqtomaf){strat <- ifelse(strat <= 0.5, strat, 1-strat)}
    nfs <- length(fs) - 1
    intervals <- findInterval(x = strat, vec = fs)
    for(i in 1:nfs){
      assign(
        paste0("l", i),
        lam(z = z[intervals == i], expected.stat = expected.stat)
      )
    }
    for(j in 1:nfs){
      PlotQQHelper(p = pval[intervals == j], color = colors[j], cex = 0.3)
    }
    legend("bottomright", legend = 
             c("Expected (null)", "Observed", 
               paste0("MAF < ", fs[2] ," [lam = ", l1, "]"), paste0(fs[2:nfs], " < MAF < ", fs[3:(nfs+1)] ," [lam = ", mget(paste0("l", 2:nfs)), "]")), 
           pch = rep(23, 6), cex = 1.1, pt.cex = 1.5, pt.bg = c("grey25", obscol, colors))
  }else{
    legend(x = "topleft", legend = c("Expected","Observed"), pch = 23, cex = 1, pt.bg = c("black", obscol), bty = "n")
  }
  if(save){dev.off()}
}
