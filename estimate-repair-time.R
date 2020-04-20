# TODO: this should be converted to a proper 2D maximum likelihood
# fit, such that the confidence region can be determined in 2D.
# The result is unlikely to change much.

set.seed(123)

# Extract calcium data with at least one hit, and order 
# all hits increasing with time
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

hit.count.function <-  function(x,t,...){
	r <- rep(0,length(x))
	for( i in seq_along(x) ){
		hit.times <- ca.hit.times[[as.character(x[i])]]
		nhits <- length(hit.times)
		last.hit <- -Inf
		for( ht in hit.times ){
			if( ht < t[i] ){
				r[i] <- r[i] + 1
			}
		}
	}
	r
}

damage.function <- function(rate){ 
	function(x,t,...){
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
}

nLL <- function( t12.minutes=30 ){
	rate <- 60*log(2)/t12.minutes
	m0 <- coxph( Surv( end.time, fate ) ~ tt(cell),
		tt=damage.function(rate))
	as.numeric(-logLik(m0))
}

# Fit the decay rate using maximum likelihood
library(stats4)
ml <- mle( nLL, start=list( t12.minutes=4 ),
	method = "L-BFGS-B", lower=c(0) )

# Plot the ML procedure
pdf("plots/estimate-repair-time.pdf", width=4, height=4)
par(font.main=1, cex.main=1)
x <- seq( 1, 240, by=5 )
y <- sapply( x, nLL )
plot(x,y, xlab="minutes", 
	xlim=c(0,180), ylab="negative log-likelihood",
	type='l', bty='l', main="estimation of recovery t1/2" )
mi <- which.max( -y )
ci <- confint(ml)
segments( ci[1], ml@min, ci[2], ml@min, col=2, lwd=2 )
points( ml@coef, ml@min, pch=19, col=2 )
dev.off()

# Compare fit of damage decay model to simple "hit counting" model

cat( "Damage model fit:\n" )
cat( "Negative log-likelihood at minimum: ", nLL( ml@coef ), "\n" )
# We need to add penalty for the extra parameter that we fitted,
# so we compute the BIC manually
bic.damage.model <- 2*log(sum(fate==1))+2*nLL( ml@coef )
cat( "BIC: ", bic.damage.model,"\n" )

cat( "\n\nSimple model fit:\n" )
mdl <- coxph( Surv(end.time,fate) ~ tt(cell), tt=hit.count.function )
bic.n.ca.hits.model <- BIC( mdl )
cat( "BIC: ", bic.n.ca.hits.model, "\n" )

cat( "\n\nBIC difference: ", bic.n.ca.hits.model-bic.damage.model,"\n" )


