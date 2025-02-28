readInLengths.fn <- function(file,headerRow=1,newNames=TRUE) {
    #make a csv file from lenght sheets
    #headerRow is the row number of the column where the data start
    #I put in simplified names. Make sure that these match
    #written by Allan Hicks, 3/15/11
    nombres <- c("ProjectCycle","Vessel","HaulID","Lat","Lon","PosType","Depth","DepthType","Sex","Length")

    xx <- read.csv(file,skip=headerRow-1)
    if(newNames) {
        names(xx) <- nombres
        cat("\nNOTE: column names have been modified from the csv file. You may want to verify that they match.\n\n")
    }
    return(xx)
}

classifyINPFC.fn <- function(x) {
    #classifies INPFC areas from latitude coordinates
    INPFClats <- c(36,40.5,43,47.5)
    INPFC <- c("CP","MT","EK","CL","VN")
    x <- findInterval(x,INPFClats)
    return(ordered(INPFC[x+1],INPFC))
}    


########################################################################################


PReg.obj <- function(x,Depth,Length,rngSplit=c(75,180),fixSplt=NULL,verbose=F) {
    #piecewise regression, assuming continuous function at breakpoint. 
    #Therefore, same intercept for both equations
    splt <- x[1]
    alpha <- x[2]
    b1 <- x[3]
    b2 <- x[4]
    nobs <- length(Depth)
    
    predY <- rep(NA,nobs)
    ind <- Depth <= splt
    predY[ind] <- alpha+b1*(Depth[ind]-splt)
    f.left <- sum((predY[ind]-Length[ind])^2)
    ind <- Depth > splt
    predY[ind] <- alpha+b2*(Depth[ind]-splt)
    f.right <- sum((predY[ind]-Length[ind])^2)

    f <- (f.left+f.right)
    if(verbose){cat(f.left,f.right,f,"\n")}
    
    if(splt<rngSplit[1] | splt>rngSplit[2]) {
        return(1e10)
    }
    return((nobs/2)*log(f/nobs))    #the likelihood function (assuming normal errors and constant variance)
}



PRegLikeProf.fn <- function(x,Depth,Length,rngSplit,step) {
  PRegLP.obj <- function(x,breakpt,Depth,Length) {
    #piecewise regression, assuming continuous function at breakpoint. 
    #Therefore, same intercept for both equations
    #Likelihood profile over splt
    splt <- breakpt
    alpha <- x[1]
    b1 <- x[2]
    b2 <- x[3]
    nobs <- length(Depth)
    predY <- rep(NA,nobs)
    ind <- Depth <= splt
    predY[ind] <- alpha+b1*(Depth[ind]-splt)
    f.left <- sum((predY[ind]-Length[ind])^2)
    ind <- Depth > splt
    predY[ind] <- alpha+b2*(Depth[ind]-splt)
    f.right <- sum((predY[ind]-Length[ind])^2)
    f <- (f.left+f.right)
    return((nobs/2)*log(f/nobs))    #the likelihood function (assuming normal errors and constant variance)
  }

    nobs <- length(Depth)
    splits <- seq(rngSplit[1],rngSplit[2],by=step)
    out <- matrix(NA,nrow=length(splits),ncol=5,dimnames=list(NULL,c("split","alpha","b1","b2","f")))
    for(i in 1:length(splits)) {
        opt <- optim(x[-1],PRegLP.obj,breakpt=splits[i],Depth=Depth,Length=Length)
        out[i,] <- c(splits[i],opt$par,opt$value)
    }
    as.data.frame(out)
}

bootstrapPReg <- function(x,Depth,Length,nBoots=1000) {
    Bout <- matrix(NA,nrow=nBoots,ncol=5,dimnames=list(NULL,c("split","alpha","b1","b2","f")))
    for(i in 1:nBoots) {
        Bind <- sample(1:length(Depth),replace=T)
        out <- optim(x,PReg.obj,Depth=Depth[Bind],Length=Length[Bind])
        Bout[i,] <- c(out$par,out$value)
        if(i%%100==0)cat("Finished",i,"simulations out of",nBoots,"\n")
    }
    as.data.frame(Bout)
}


PReg.fn <- function(Depth,Length,rngSplit=c(55,183),step=1,const=100) {
    nobs <- length(Depth)
    splits <- seq(rngSplit[1],rngSplit[2],by=step)
    out <- matrix(NA,nrow=length(splits),ncol=6,dimnames=list(NULL,c("split","a1","b1","a2","b2","f")))
    for(i in 1:length(splits)) { 
        splt <- splits[i]
        predY <- rep(NA,nobs)
        ind <- Depth <= splt
        tmp1.lm <- lm(Length[ind]~Depth[ind])
        f <- sum((tmp1.lm$residuals/const)^2)
        ind <- Depth > splt
        tmp2.lm <- lm(Length[ind]~Depth[ind])
        f <- f+sum((tmp2.lm$residuals/const)^2)
        out[i,] <- c(splt,tmp1.lm$coefficients,tmp2.lm$coefficients,(nobs/2)*log(f/nobs))
    }
    as.data.frame(out)
}



#######################################################################################################
#uses tree regression to determine strata (see Francis FAR 2002/9 and Hicks FAR 2002/?)
treeStrata.fn <- function(obs,responseVar,splitVars,minNumInStrata=10,minPercentS0=1,cohort=NA,beta=NULL) {
    #Uses a dataframe with the variables specified
    #   obs: the dataframe of data
    #   responseVar: the name of the response variable (column name in obs)
    #minNumStrata is the minimum number of observations per strata
    #minSS is the minimum additional percent of the unstratified sum of squares to explain, i.e. 100*((S[k]-S[k+1])/S[0])

    #this function will print out a tree, cut back to satisfy the criteria entered
    #it will return a tree object with the full tree (untrimmed)
    #Each row of the printed tree contains the following:
    # variable split on, number in split, deviance, mean length, (percent additional variation explained by split below)
    #if the message "CHECK THE ADDITIONAL VARAINCE EXPLAINED CRITERIA BY HAND" appears, you should look at the full tree
      #and trim it back accordingly because it was a bit to complicated for this simple program
      #one way to do this is to look at the tree with your minPercentS0 criteria and take a look at the tree with
      #the minPercentS0 criteria set to zero (full tree) and see that no important splits were removed
      #doing it this way will print out the percent additional variation explained for the full tree
    #I suggest that when a strata says less than ("<"), it should be interpreted as less than or equal to
    
    obs <- data.frame(response=obs[,responseVar],obs[,splitVars])
    stratTree <- tree(response~.,data=obs,mincut=minNumInStrata)
    fullTree <- stratTree
    
    tempfr <- stratTree$frame
    tempfr$percAddDev <- rep(NA,nrow(tempfr)) 
    S0 <- tempfr[1,"dev"]
    for(i in 1:nrow(tempfr)) {
        if(tempfr[i,"var"] != "<leaf>") {
            nextNode <- as.numeric(row.names(tempfr)[i+1])
            nodes <- as.character(c(nextNode,nextNode+1))
            tempfr$percAddDev[i] <- 100*(tempfr[i,"dev"]-sum(tempfr[nodes,"dev"]))/S0
        }
    }
    noSplits <- tempfr$percAddDev<minPercentS0
    noSplits[is.na(noSplits)] <- F
    noSplits <- row.names(tempfr[noSplits,])
    if(length(noSplits)>0) {
        for(ii in length(noSplits):1) {
            temp.nodes <- as.character(c(as.numeric(noSplits[ii])*2,as.numeric(noSplits[ii])*2+1))
            if(tempfr[temp.nodes[1],"var"] == "<leaf>" && tempfr[temp.nodes[2],"var"] == "<leaf>") {
                rows <- (1:nrow(tempfr))[row.names(tempfr)==temp.nodes[1]]
                tempfr[rows-1,"var"] <- "<leaf>"
                tempfr[rows-1,"percAddDev"] <- NA
                tempfr[rows-1,"splits"] <- ""
                rows <- c(rows,rows+1)
                tempfr <- tempfr[-rows,]
            }
            else { #don't delete because there is a lower split
                print("CHECK THE ADDITIONAL VARAINCE EXPLAINED CRITERIA BY HAND")
            }
        }
    }
    else{
        print("No splits were removed")
    }
    if(nrow(tempfr)>1) {
        stratTree$frame <- tempfr
        print.tree.ach(stratTree)
    }
    else {
        print("There are no suitable strata")
    }
    fullTree
}




print.tree.ach <- function(x, pretty = 0, spaces = 2, digits = .Options$digits - 3, ...)
{
#spaces is used for indentation of tree levels
    if(!inherits(x, "tree")) stop("Not legitimate tree")
    if(is.prob <- !is.null(ylevels <- attr(x, "ylevels")))
        cat("node), split, n, deviance, yval, (yprob)\n")
    else cat("node), split, n, deviance, yval (percAddVar of split below)\n")
    cat("      * denotes terminal node\n\n")
    frame <- x$frame
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    indent <- paste(rep(" ", spaces * 32), collapse = "")   #32 is the maximal depth
    if(length(node) > 1) {
        indent <- substring(indent, 1, spaces * seq(depth))
        indent <- paste(c("", indent[depth]), format(node), ")", sep = "")
    }
    else indent <- paste(format(node), ")", sep = "")
    if(is.prob) {
        yval <- paste(as.character(frame$yval), " (", sep = "")
        yprob <- frame$yprob
        for(i in seq(ylevels))
            yval <- paste(yval, format(signif(yprob[, i], digits = digits)))
        yval <- paste(yval, ")")
    }
    else yval <- format(signif(frame$yval, digits = digits))
    term <- rep(" ", length(depth))
    term[frame$var == "<leaf>"] <- "*"
    z <- labels.tree(x, pretty = pretty)
    if(inherits(x, "bonzai"))
        n <- format(signif(frame$n, digits = digits))
    else n <- frame$n
    z <- paste(indent, z, n, format(signif(frame$dev, digits = digits)), yval, paste("(",round(frame$percAddDev,1),")",sep=""), term)
    cat(z, sep = "\n")
    invisible(x)
}





tree.plot.ach <- function (x, y = NULL, type = c("proportional", "uniform"),label=T, ...) 
{
    if (inherits(x, "singlenode")) 
           stop("cannot plot singlenode tree")
    if (!inherits(x, "tree")) 
        stop("not legitimate tree")
    type <- match.arg(type)
    uniform <- type == "uniform"
    dev <- dev.cur()
    if (dev == 1) 
        dev <- 2
    assign(paste(".Tree.unif", dev, sep = "."), uniform, envir = .GlobalEnv)
    tmp <- treeplRotate(treeco(x), node = as.numeric(row.names(x$frame)),...)
    if(is.logical(label)) {
        if(!label)   return(tmp)
    }

    #otherwise plot the labels
    xy <- data.frame(tmp,x$frame)
    xy$depth <- tree.depth(as.numeric(row.names(x$frame)))
    splits <- xy$splits[xy$splits[,1]!="",]
    xyt <- xy[xy$var!="<leaf>",names(xy)!="splits"]
    xyt <- data.frame(rbind(xyt,xyt),splits=as.vector(splits))
    if(!is.logical(label)) {
        xyt$label <- as.character(label)
    }    
    else {
        xyt$label <- paste(xyt$var,xyt$splits)
        ind <- xyt$var=="Year"
        xyt$label[ind] <- paste("Year",substring(xyt$splits[ind],1,1),floor(as.numeric(substring(xyt$splits[ind],2))),sep="")
        ###can do some more default stuff
    }
    xyt$z <- rep(NA,nrow(xyt))
    #xyt <- xyt[order(xyt$x,xyt$y),]
    #xyt$z[
    print(xyt)
    charht <- par("cxy")[2]
    text(xyt$x,xyt$y+0.5*charht,xyt$label,adj=0)
    
    

}
#tree.plot.ach(treeB)


treepl <- function (xy, node, erase = FALSE, ...) 
{
    x <- xy$x
    y <- xy$y
    parent <- match((node%/%2), node)
    sibling <- match(ifelse(node%%2, node - 1, node + 1), node)
    xx <- rbind(x, x, x[sibling], x[sibling], NA)
    yy <- rbind(y, y[parent], y[parent], y[sibling], NA)
    if (any(erase)) {
        lines(c(xx[, erase]), c(yy[, erase]), col = par("bg"))
        return(x = x[!erase], y = y[!erase])
    }
    plot(range(x), range(y), type = "n", axes = FALSE, xlab = "", 
        ylab = "")
    text(x[1], y[1], "|", ...)
    lines(c(xx[, -1]), c(yy[, -1])-0.5*par("cxy")[2], ...)
    list(x = x, y = y)
}

treeplRotate <- function (xy, node, erase = FALSE, ...) 
{
    x <- xy$x
    y <- -1*xy$y
    parent <- match((node%/%2), node)
    sibling <- match(ifelse(node%%2, node - 1, node + 1), node)
    xx <- rbind(x, x, x[sibling], x[sibling], NA)
    yy <- rbind(y, y[parent], y[parent], y[sibling], NA)
    if (any(erase)) {
        lines(c(xx[, erase]), c(yy[, erase]), col = par("bg"))
        return(x = x[!erase], y = y[!erase])
    }
    plot(range(y), range(x), type = "n", axes = FALSE, xlab = "", 
        ylab = "")
    text(y[1], x[1], "-", ...)
    lines(c(yy[, -1]), c(xx[, -1]), ...)
    list(x = y, y = x)
}

treeco <- function (tree, uniform = paste(".Tree.unif", dev.cur(), sep = ".")) 
{
    frame <- tree$frame
    node <- as.numeric(row.names(frame))
    depth <- tree.depth(node)
    x <- -depth
    if (exists(uniform)) 
        uniform <- get(uniform)
    else uniform <- 0
    if (uniform) 
        y <- x
    else {
        y <- dev <- frame$dev
        depth <- split(seq(node), depth)
        parent <- match(node%/%2, node)
        sibling <- match(ifelse(node%%2, node - 1, node + 1), 
            node)
        for (i in depth[-1]) y[i] <- y[parent[i]] - dev[parent[i]] + 
            dev[i] + dev[sibling[i]]
    }
    depth <- -x
    leaves <- frame$var == "<leaf>"
    x[leaves] <- seq(sum(leaves))
    depth <- split(seq(node)[!leaves], depth[!leaves])
    left.child <- match(node * 2, node)
    right.child <- match(node * 2 + 1, node)
    for (i in rev(depth)) x[i] <- 0.5 * (x[left.child[i]] + x[right.child[i]])
    list(x = x, y = y)
}

tree.depth <- function (nodes) 
{
    depth <- floor(log(nodes, base = 2) + 1e-07)
    as.vector(depth - min(depth))
}

text.treeRotate <- function (x, splits = TRUE, label = "yval", all = FALSE, pretty = NULL, 
    digits = getOption("digits") - 3, adj = par("adj"), xpd = TRUE, 
    ...) 
{
    oldxpd <- par(xpd = xpd)
    on.exit(par(oldxpd))
    if (inherits(x, "singlenode")) 
        stop("cannot plot singlenode tree")
    if (!inherits(x, "tree")) 
        stop("not legitimate tree")
    frame <- x$frame
    column <- names(frame)
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
        column <- c(column, ylevels)
    if (!is.null(label) && is.na(match(label, column))) 
        stop("label must be a column label of the frame component of the tree")
    charht <- par("cxy")[2]
    if (!is.null(srt <- list(...)$srt) && srt == 90) {
        if (missing(adj)) 
            adj <- 0
        ladj <- 1 - adj
    }
    else ladj <- adj
    xy <- treeco(x)
    yy <- xy$x
    xx <- -1*xy$y

    if (splits) {
        node <- as.numeric(row.names(frame))
        left.child <- match(2 * node, node)
        rows <- labels.tree(x, pretty = pretty)[left.child]
        ind <- !is.na(rows)
        text(xx[ind], yy[ind] + 0.5 * charht, rows[ind], 
            adj = adj, ...)
    }
    if (!is.null(label)) {
        leaves <- if (all) 
            rep(TRUE, nrow(frame))
        else frame$var == "<leaf>"
        if (label == "yval" & !is.null(ylevels)) 
            stat <- as.character(frame$yval[leaves])
        else if (!is.null(ylevels) && !is.na(lev <- match(label, 
            ylevels))) 
            stat <- format(signif(frame$yprob[leaves, lev], digits = digits))
        else stat <- format(signif(frame[leaves, label], digits = digits))
        if (!is.null(dim(stat)) && dim(stat)[2] > 1) {
            if (length(dimnames(stat)[[2]])) 
                stat[1, ] <- paste(sep = ":", dimnames(stat)[[2]], 
                  stat[1, ])
            stat <- do.call("paste", c(list(sep = "\n"), split(stat, 
                col(stat))))
        }
        text(xx[leaves], yy[leaves] - 0.5 * charht, labels = stat, 
            adj = ladj, ...)
    }
    invisible()
}


#text.treeRotate(treeB)







text.tree <- function (x, splits = TRUE, label = "yval", all = FALSE, pretty = NULL, 
    digits = getOption("digits") - 3, adj = par("adj"), xpd = TRUE, 
    ...) 
{
    oldxpd <- par(xpd = xpd)
    on.exit(par(oldxpd))
    if (inherits(x, "singlenode")) 
        stop("cannot plot singlenode tree")
    if (!inherits(x, "tree")) 
        stop("not legitimate tree")
    frame <- x$frame
    column <- names(frame)
    if (!is.null(ylevels <- attr(x, "ylevels"))) 
        column <- c(column, ylevels)
    if (!is.null(label) && is.na(match(label, column))) 
        stop("label must be a column label of the frame component of the tree")
    charht <- par("cxy")[2]
    if (!is.null(srt <- list(...)$srt) && srt == 90) {
        if (missing(adj)) 
            adj <- 0
        ladj <- 1 - adj
    }
    else ladj <- adj
    xy <- treeco(x)
    if (splits) {
        node <- as.numeric(row.names(frame))
        left.child <- match(2 * node, node)
        rows <- labels.tree(x, pretty = pretty)[left.child]
        ind <- !is.na(rows)
        text(xy$x[ind], xy$y[ind] + 0.5 * charht, rows[ind], 
            adj = adj, ...)
    }
    if (!is.null(label)) {
        leaves <- if (all) 
            rep(TRUE, nrow(frame))
        else frame$var == "<leaf>"
        if (label == "yval" & !is.null(ylevels)) 
            stat <- as.character(frame$yval[leaves])
        else if (!is.null(ylevels) && !is.na(lev <- match(label, 
            ylevels))) 
            stat <- format(signif(frame$yprob[leaves, lev], digits = digits))
        else stat <- format(signif(frame[leaves, label], digits = digits))
        if (!is.null(dim(stat)) && dim(stat)[2] > 1) {
            if (length(dimnames(stat)[[2]])) 
                stat[1, ] <- paste(sep = ":", dimnames(stat)[[2]], 
                  stat[1, ])
            stat <- do.call("paste", c(list(sep = "\n"), split(stat, 
                col(stat))))
        }
        text(xy$x[leaves], xy$y[leaves] - 0.5 * charht, labels = stat, 
            adj = ladj, ...)
    }
    invisible()
}


#From old libaries of tree
labels.tree <- function (object, pretty = TRUE, collapse = TRUE, ...) 
{
    if (!inherits(object, "tree")) 
        stop("not legitimate tree")
    frame <- object$frame
    xlevels <- attr(object, "xlevels")
    var <- as.character(frame$var)
    splits <- matrix(sub("^>", " > ", sub("^<", " < ", frame$splits)), 
        , 2)
    if (!is.null(pretty)) {
        if (pretty) 
            xlevels <- lapply(xlevels, abbreviate, minlength = pretty)
        for (i in grep("^:", splits[, 1])) for (j in 1L:2L) {
            sh <- splits[i, j]
            nc <- nchar(sh)
            sh <- substring(sh, 2L:nc, 2L:nc)
            xl <- xlevels[[var[i]]][match(sh, letters)]
            splits[i, j] <- paste(": ", paste(as.vector(xl), 
                collapse = ","), sep = "")
        }
    }
    if (!collapse) 
        return(array(paste(var, splits, sep = ""), dim(splits)))
    node <- as.numeric(row.names(frame))
    parent <- match((node%/%2), node)
    odd <- as.logical(node%%2)
    node[odd] <- paste(var[parent[odd]], splits[parent[odd], 
        2L], sep = "")
    node[!odd] <- paste(var[parent[!odd]], splits[parent[!odd], 
        1L], sep = "")
    node[1L] <- "root"
    node
}
