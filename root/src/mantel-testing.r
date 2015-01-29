library(igraph)
library(maps)
library(vegan)
library(RColorBrewer)
library(fields)

## Get case data

dataDir <- '/root/data'

tmpf <- function(){
  fn <- file.path(dataDir, 'PEDvweeklyreport-state-ts-01-08-14.csv')
  ret <- read.csv(fn)
  ret[1:9, c('CA', 'MD', 'NE', 'WY')] <- 0
  key <- order(colnames(ret)[-c(1:2)])
  ret <- cbind(ret[, c(1:2)], ret[, -c(1:2)][, key])
  ret$Unk <- NULL
  target <- structure(list(week = structure(c(20L, 29L, 32L, 8L, 12L),
.Label = c("10/13/2013", "10/20/2013", "10/27/2013", "10/6/2013",
"11/10/2013", "11/17/2013", "11/24/2013", "11/3/2013", "12/1/2013",
"12/15/2013", "12/22/2013", "12/29/2013", "12/8/2013", "4/15/2013",
"4/22/2013", "4/29/2013", "5/13/2013", "5/20/2013", "5/27/2013",
"5/6/2013", "6/10/2013", "6/16/2013", "6/23/2013", "6/3/2013",
"6/30/2013", "7/14/2013", "7/21/2013", "7/28/2013", "7/7/2013",
"8/11/2013", "8/18/2013", "8/25/2013", "8/4/2013", "9/1/2013",
"9/15/2013", "9/22/2013", "9/29/2013", "9/8/2013"), class = "factor"),
totalNumberSwineAccessions = c(17L, 34L, 26L, 90L, 134L), CA = c(0, 0,
0, 0, 1), CO = c(1L, 1L, 0L, 0L, 0L), IA = c(8L, 6L, 2L, 38L, 54L), IL
= c(0L, 1L, 0L, 1L, 14L), IN = c(3L, 2L, 1L, 1L, 3L), KS = c(0L, 4L,
3L, 6L, 4L), KY = c(0L, 0L, 0L, 0L, 0L), MD = c(0, 0, 0, 0, 0), MI =
c(0L, 0L, 0L, 2L, 0L), MN = c(1L, 2L, 2L, 7L, 20L), MO = c(0L, 0L, 0L,
2L, 4L), NC = c(0L, 3L, 4L, 14L, 18L), NE = c(0, 0, 0, 0, 2), NY =
c(0L, 0L, 0L, 0L, 0L), OH = c(0L, 2L, 1L, 5L, 5L), OK = c(0L, 11L,
10L, 10L, 2L), PA = c(1L, 0L, 3L, 1L, 0L), SD = c(0L, 0L, 0L, 0L, 3L),
TN = c(0L, 0L, 0L, 1L, 1L), TX = c(0L, 0L, 0L, 1L, 1L), WI = c(0L, 0L,
0L, 1L, 0L ), WY = c(0, 0, 0, 0, 1)), .Names = c("week",
"totalNumberSwineAccessions", "CA", "CO", "IA", "IL", "IN", "KS",
"KY", "MD", "MI", "MN", "MO", "NC", "NE", "NY", "OH", "OK", "PA",
"SD", "TN", "TX", "WI", "WY" ), row.names = c(4L, 13L, 20L, 30L, 38L),
class = "data.frame")
  stopifnot(identical(target, ret[c(4,13,20,30,38),]))
  ret
}

caseData <- tmpf()

unwanted <- c('week', 'totalNumberSwineAccessions', 'Unk')
ind <- which(!colnames(caseData) %in% unwanted)
observed <- caseData[, ind]

## Get cross correlations

n <- ncol(observed)
CC <- matrix(nrow=n, ncol=n)

getcc <- function(x,y, lag=1){
    foo <- ccf(x, y, plot=FALSE)
    ind <- which(foo$lag == lag)
    foo$acf[ind]
}

for(i in seq_len(n)){
    for(j in seq_len(n)){
        ## CC[i,j] will be high if deviations from the mean in series i
        ## are shifted to the left of similar deviations to the mean in series j
        ## i.e. i's deviations are indicative of j's future deviations
        CC[i,j] <- getcc(observed[,j], observed[,i])
    }
}

colnames(CC) <- rownames(CC) <- colnames(observed)

check_cc_directionality <- function(){
    t <- 1:100 * .25
    x <- sin(t)
    ## y is shifted to the right
    y <- sin(t-1)

    xylag1 <- getcc(x, y)
    yxlag1 <- getcc(y,x)
    main <- paste('getcc(x, y) ==', round(xylag1,2),
                  '; getcc(y,x) == ', round(yxlag1,2))
    if(xylag1 > yxlag1){
        conc <- 'left arg delayed by lag'
    } else{
        conc <- 'right arg delayed by lag'
    }
    plot(x~t, type='l', ylab='f(t)', main=main, sub=conc)
    lines(y~t, col=2)
    legend('topright', col=1:2, legend=c('x', 'y'), lty=1)
}

png('direction-check.png')
check_cc_directionality()
dev.off()

## Get shared-border neighborhoods

nbEdgelist <- read.csv(file.path(dataDir, 'state_neighbors_fips.txt'), header=FALSE)
data(state.fips)
g <- graph.data.frame(nbEdgelist, directed=FALSE)
g <- simplify(g)
key <- match(V(g)$name, state.fips$fips)
abb <- state.fips$abb[key]
V(g)$name <- as.character(abb)
vids <- which(V(g)$name %in% colnames(caseData))
g2 <- induced.subgraph(g, vids)
nms <- colnames(observed)
nhood <- as.matrix(g2[nms,nms])

## Get shipment flows

ep <- read.csv(file.path(dataDir, 'edge-proportionalities.csv'))
key <- match(colnames(observed), colnames(ep))
ep <- ep[key, key]
rownames(ep) <- colnames(ep)

### assertions based on inspection of original .xls file
stopifnot(ep['CA', 'IL'] == 9415)
stopifnot(ep['MI', 'KS'] == 6)
stopifnot(ep['IL', 'IA'] == 1424813)
stopifnot(ep['KY', 'IA'] == 17658)
stopifnot(ep['MO', 'IA'] == 2389932)
stopifnot(ep['OK', 'CA'] == 16762)
stopifnot(ep['MO', 'CA'] == 645)

epl <- log10(ep +1)
epl <- data.matrix(epl)
## We discard the within-state flows because they do not enter into
## the analysis. Thus it is mostly a matter of how the plots will
## look, and because the original spreadsheet had no data on
## within-state flows (they were added for some other analyses based
## on an assumption of 90% of total flow), there's no good reason to
## include them in the plot.
diag(epl) <- NA

key <- match(colnames(ep), colnames(CC))
CC <- CC[key,key]



## Get great circle distance

key <- match(colnames(CC), state.abb)

cx <- state.center$x[key]
cy <- state.center$y[key]

n <- ncol(CC)
centerDists <- matrix(nrow=n, ncol=n)

# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# Haversine formula (hf)
# source: http://www.r-bloggers.com/great-circle-distance-calculations-in-r/
gcd.hf <- function(long1, lat1, long2, lat2) {
      R <- 6371 # Earth mean radius [km]
        delta.long <- (long2 - long1)
        delta.lat <- (lat2 - lat1)
        a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
        c <- 2 * asin(min(1,sqrt(a)))
        d = R * c
        return(d) # Distance in km
  }

deg2rad <- function(deg) return(deg*pi/180)


getCenterDist <- function(state1,state2){
    long1 <- deg2rad(cx[state1])
    long2 <- deg2rad(cx[state2])
    lat1 <- deg2rad(cy[state1])
    lat2 <- deg2rad(cy[state2])
    gcd.hf(long1, lat1, long2, lat2)
}

for(i in seq_len(n)){
    for(j in seq_len(n)){
        centerDists[i,j] <- getCenterDist(i, j)
    }
}

colnames(centerDists) <- rownames(centerDists) <- colnames(CC)

## hypothesis tesing

doTest <- function(M1, M2, symmetrize=FALSE, ...){
    getMat <- function(x) switch(x, 'shipment'=epl, 'cor'=CC,
                     'gcd'=-centerDists, 'sharedBord'=nhood)
    x <- getMat(M1)
    y <- getMat(M2)
    if(symmetrize){
        x <- x + t(x)
        y <- y + t(y)
    }
    mantel(x, y, ...)    
}

mats <- c('shipment', 'cor', 'gcd', 'sharedBord')
methods <- c('spearman', 'pearson')
symmetrize <- c(TRUE, FALSE)
des <- expand.grid(M1=mats, M2=mats, method=methods, symmetrize=symmetrize, stringsAsFactors=FALSE)
des <- des[des$M1 != des$M2, ]
des$permutations <- 10000

res <- list()
for(i in seq_len(nrow(des))){
    print(i)
    res[[i]] <- do.call(doTest, as.list(des[i,]))
}

des$r <- sapply(res, '[[', 'statistic')
des$pValues <- sapply(res, '[[', 'signif')

sink('mantel-table.txt')
pander::pander(des)
sink()

## plotting

### diagnostics

df <- data.frame(epl=unclass(as.dist(epl)), CC=unclass(as.dist(CC)))
m <- lm(CC~epl, data=df)
pdf('diagnostics.pdf')
plot(m)
dev.off()

### correlelogram                   

#mcor <- mantel.correlog(D.eco=CC, D.geo=epl, r.type='pearson')


### heat map

image.plot.ebo <- function (..., add = FALSE, nlevel = 64, horizontal = FALSE,
                            legend.shrink = 0.9, legend.width = 1.2,
                            legend.mar = ifelse(horizontal, 3.1, 5.1),
                            legend.lab = NULL, legend.line = 2,
                            graphics.reset = FALSE, bigplot = NULL,
                            smallplot = NULL, legend.only = FALSE,
                            col = tim.colors(nlevel), lab.breaks = NULL,
                            axis.args = NULL, legend.args = NULL,
                            midpoint = FALSE, border = NA, lwd = 1,
                            panelLab=NULL) {
    old.par <- par(no.readonly = TRUE)
    info <- imageplot.info(...)
    if (add) {
        big.plot <- old.par$plt
    }
    if (legend.only) {
        graphics.reset <- TRUE
    }
    if (is.null(legend.mar)) {
        legend.mar <- ifelse(horizontal, 3.1, 5.1)
    }
    temp <- imageplot.setup(add = add, legend.shrink = legend.shrink, 
                            legend.width = legend.width, legend.mar = legend.mar, 
                            horizontal = horizontal, bigplot = bigplot, smallplot = smallplot)
    smallplot <- temp$smallplot
    bigplot <- temp$bigplot
    if (!legend.only) {
        if (!add) {
            par(plt = bigplot)
        }
        if (!info$poly.grid) {
            image(..., add = add, col = col, axes=FALSE)
            addAxis1 <- function(x, y, z,...){
                labs <- paste(names(y), c('', '      '))
                axis(1, at=y, labels = labs, las = 2, line = -0.5, tick = 0, 
        cex.axis = 0.8)
            }
            addAxis1(...)
            addAxis2 <- function(x, y, z,...){
                labs <- paste(names(x), c('', '      '))
                axis(2, at=x, labels = labs, las = 1, line = -0.5, tick = 0, 
        cex.axis = 0.8)
            }
            addAxis2(...)
            mtext(panelLab, side=2, las=2, at=par()$usr[2], line=par()$mar[2] - 1, cex=1.5)
        }
        else {
            poly.image(..., add = add, col = col, midpoint = midpoint, 
                       border = border, lwd.poly = lwd)
        }
        big.par <- par(no.readonly = TRUE)
    }
    if ((smallplot[2] < smallplot[1]) | (smallplot[4] < smallplot[3])) {
        par(old.par)
        stop("plot region too small to add legend\n")
    }
    ix <- 1
    minz <- info$zlim[1]
    maxz <- info$zlim[2]
    binwidth <- (maxz - minz)/nlevel
    midpoints <- seq(minz + binwidth/2, maxz - binwidth/2, by = binwidth)
    iy <- midpoints
    iz <- matrix(iy, nrow = 1, ncol = length(iy))
    breaks <- list(...)$breaks
    par(new = TRUE, pty = "m", plt = smallplot, err = -1)
    if (!is.null(breaks) & !is.null(lab.breaks)) {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2), 
                            at = breaks, labels = lab.breaks), axis.args)
    }
    else {
        axis.args <- c(list(side = ifelse(horizontal, 1, 4), 
                            mgp = c(3, 1, 0), las = ifelse(horizontal, 0, 2)), 
                       axis.args)
    }
    if (!horizontal) {
        if (is.null(breaks)) {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col)
        }
        else {
            image(ix, iy, iz, xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col, breaks = breaks)
        }
    }
    else {
        if (is.null(breaks)) {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col)
        }
        else {
            image(iy, ix, t(iz), xaxt = "n", yaxt = "n", xlab = "", 
                  ylab = "", col = col, breaks = breaks)
        }
    }
    do.call("axis", axis.args)
    box()
    if (!is.null(legend.lab)) {
        legend.args <- list(text = legend.lab, side = ifelse(horizontal, 
                                                   1, 4), line = legend.line)
    }
    if (!is.null(legend.args)) {
        do.call(mtext, legend.args)
    }
    mfg.save <- par()$mfg
    if (graphics.reset | add) {
        par(old.par)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
    else {
        par(big.par)
        par(plt = big.par$plt, xpd = FALSE)
        par(mfg = mfg.save, new = FALSE)
        invisible()
    }
}

pdf('cors-heatmap.pdf')
orderings <- heatmap(data.matrix(epl), scale='none')
dev.off()

makeImagePlot <- function(M, orderings, ...){
    x <- data.matrix(M)[orderings$rowInd, orderings$colInd]
    x <- t(x)
    col <- brewer.pal(9, 'Blues')
    nc <- ncol(x)
    nr <- nrow(x)
    xi <- 1:nr
    yi <- 1:nc
    names(xi) <- colnames(x)
    names(yi) <- rownames(x)
    image.plot.ebo(x=xi, y=yi, z=x, horizontal=FALSE,col=col, graphics.reset=TRUE, ...)
}


pdf('matrices.pdf', width=8.7/2.54, height=12/2.54)
layout(matrix(1:2, ncol=1))
par(mar=c(4,4.5,.5,.5))
makeImagePlot(M=data.matrix(epl), orderings=orderings, xlab='Destination', ylab='Source',
              legend.args=list(text=expression(paste(Log[10], '(swine shipped)')),
                  line=2.9, side=4), panelLab='A')
par(mar=c(4,4.5,1,.5))
CCnoDiag <- CC
diag(CCnoDiag) <- NA
makeImagePlot(M=data.matrix(CCnoDiag), orderings=orderings, xlab='Leading state', ylab='Lagging state', 
              legend.args=list(text='Cross correlation', line=2.9, side=4), panelLab='B')
dev.off()

##

library(ggplot2)
library(reshape2)
library(scales)
library(plyr)

mts <- melt(caseData, id='week')
keepers <- c('MN', 'KS',
             'IL', 'OK',
             'IA', 'NC')
test <- mts$variable %in% keepers
mts <- mts[test, ]
mts$variable <- factor(mts$variable, levels=keepers)
mts$x <- as.Date(as.character(mts$week), format='%m/%d/%Y')
labdf <- ddply(mts, 'variable', summarize, minx=x[1], maxy=max(value))
labdf$minx[labdf$variable == 'OK'] <- as.Date('2013-11-01')

tmpf <- function(df){
    g <- ggplot(df, aes(x=x, y=value, group=variable))
    g <- g + geom_hline(yintercept=0, size=0.5, col='grey')
    g <- g + geom_step(direction="vh")
    g <- g + facet_wrap(~variable, ncol=2, scales='free_y')
    g <- g + scale_x_date()
    g <- g + scale_y_discrete(breaks=pretty_breaks(n=2))
    g <- g + xlab('Date') + ylab('Cases')
    g <- g + theme_classic()
    g <- g + theme(strip.background = element_blank(),
                   strip.text.x = element_blank())
    g <- g + theme(plot.margin=unit(c(0,2,0,0),"mm"))
    g <- g + geom_text(data=labdf, hjust=0, vjust=1,
                       aes(x=minx, y=maxy, label=variable))
    g
}

g <- tmpf(mts)

ggsave('ts.pdf', width=8.6/2.54, height=6/2.54) 
