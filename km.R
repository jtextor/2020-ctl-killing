

d <- read.table("data/calcium-trace.txt", header=TRUE)[,-2]
d$event.type[d$event.type==-1] <- 0

# Extract only those cells that have had at least 1 Ca2+ event
nr <- c( by(d, d$cell, nrow) )
d <- d[d$cell %in% names(which(nr>1)),]


plt <- function(d){
	library( survival )
	lag.time <- c( by( d$time, d$cell, function(x) diff(tail(x,2)) ) )
	fate <- c( by(d$event.type, d$cell, function(x) tail(x,1)) )
	prior.hits <- c( by( d, d$cell, nrow ) )-2
	plot( survfit( Surv(lag.time[prior.hits==0], fate[prior.hits==0]) ~ 1 ),
		conf.int=FALSE, xlab="time after last Ca2+ hit (h)", ylab="survival" )
	lines( survfit( Surv(lag.time[prior.hits==1], fate[prior.hits==1]) ~ 1 ),
		conf.int=FALSE, col=2 )
	lines( survfit( Surv(lag.time[prior.hits==2], fate[prior.hits==2]) ~ 1 ),
		conf.int=FALSE, col=3 )
	lines( survfit( Surv(lag.time[prior.hits>2], fate[prior.hits>2]) ~ 1 ),
		conf.int=FALSE, col=4 )
}



pdf("plots/km.pdf", width=4, height=3 )
par(mar=c(4,4,.2,.2),font.main=1, cex.main=1, bty='l')
plt(d)
dev.off()



