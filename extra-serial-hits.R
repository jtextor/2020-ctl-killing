

d <- read.table("data/calcium-trace.txt", header=TRUE)[,-2]
d$event.type[d$event.type==-1] <- 0

# Extract only those cells that have had at least 1 Ca2+ event
nr <- c( by(d, d$cell, nrow) )
d <- d[d$cell %in% names(which(nr>1)),]

times <- c()
fates <- c()
the.min.events <- c()

recovery.time <- 56.73069 / 60


for( min.events in c(1,2,3) ){

for( i in head( unique( d$cell ), Inf ) ){
	ds <- d[d$cell==i,]
	n.events <- rep( 0, nrow(ds) )
	for( j in seq_len( nrow(ds)-1 ) ){
		n.events[j] <- sum( head(ds$time,j) >= ds$time[j] - recovery.time )
	}
	ds <- cbind( ds, n.events )
	if( any( n.events >= min.events ) ){
		first.stretch <- which( n.events >= min.events )[1]
		ds <- ds[c(first.stretch, nrow(ds)),]
		ds$time <- ds$time - ds$time[1]
		print( ds )
		times <- c( times, ds$time[2] )
		fates <- c( fates, ds$event.type[2] )
		the.min.events <- c( the.min.events, min.events ) 
	}
}

}

pdf("plots/extra-serial-hits.pdf", width=6, height=4 )
par( mar=c(4,4,1,1) )
library( survival )
fit <- survfit( Surv(times,fates) ~ the.min.events )
plot( fit , conf.int=F, mark=T, pch=19, bty="l", yscale=100,
	xlab="time (h)", ylab="survival (%)", col=1:3, xlim=c(0,6) )
# abline( v = 2 )
legend( "topright", c("1 hit","2 hits","3 hits"), lty=1, col=1:3, bty="n" )
dev.off()

