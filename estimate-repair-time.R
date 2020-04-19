# TODO: this should be converted to a proper 2D maximum likelihood
# fit, such that the confidence region can be determined in 2D.
# The result is unlikely to change much.

set.seed(123)

d <- read.table("data/calcium-trace.txt", header=TRUE)[,-2]
colnames(d) <- c("cell","time","event.type","ctl")
d$event.type[d$event.type==-1] <- 0
d <- d[order(d$cell,d$time),]

n.ca.hits <- c(by( d, d$cell, function(x)nrow(x)-1 ))
cells <- unique(d$cell)

final.event <- which( d$event.type <= 1 )

end.time <- d$time[final.event]
fate <- d$event.type[final.event]
cell <- d$cell[final.event]

library( survival )

ca.hit.times <- lapply( split( d$time, d$cell ), function(x) c(head(x,-1),Inf) )

nLL <- function( t12=30 ){
	rate <- 60*log(2)/t12
	m0 <- coxph( Surv( rep(0,length(end.time)), end.time, fate ) ~ tt(cell),
		tt=function(x,t,...){
			r <- rep(0,length(x))
			for( i in seq_along(x) ){
				hit.times <- ca.hit.times[[as.character(x[i])]]
				nhits <- length(hit.times)
				last.hit <- -Inf
				for( ht in hit.times ){
					if( ht < t[i] ){
						r[i] <- r[i] + dexp( t[i]-ht, rate )
					}
				}
			}
			r
		}
	)
	as.numeric(-logLik(m0))
}

x <- seq( 1,240, by=5)
y <- sapply( x, nLL )


library(stats4)

ml <- mle( nLL, start=list( t12=4 ),
	method = "L-BFGS-B", lower=c(0) )

pdf("plots/estimate-repair-time.pdf", width=4, height=4)

par(font.main=1, cex.main=1)

plot(x,y, xlab="minutes", 
	xlim=c(0,180), ylab="negative log-likelihood",
	type='l', bty='l', main="estimation of recovery t1/2" )

mi <- which.max( -y )
ci <- confint(ml)
segments( ci[1], ml@min, ci[2], ml@min, col=2, lwd=2 )
points( ml@coef, ml@min, pch=19, col=2 )

dev.off()

