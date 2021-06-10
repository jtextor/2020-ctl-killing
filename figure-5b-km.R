library( rms )

load.data.raw <- function(){
	# print("trying to load raw  data")
	d <- read.table("data/calcium-trace.txt", head=TRUE)

	d$event.type[d$event.type==-1] <- 0

	# Extract only those cells that have had at least 1 Ca2+ event
	nr <- c( by(d, d$cell, nrow) )
	d <- d[d$cell %in% names(which(nr>1)),]

	#calc.events <- by(d$event.type, d$cell, function(d) sum(d==2) )
	# our model assumes that death does not happen spontaneously
	#d <- d[d$cell %in% names(which(calc.events>0)),]

	oi <- unique(d$cell)
	d$cell <- sapply( d$cell, function(x) which(oi==x))
	d[order(d$cell,d$time),]
}

# explain death time with any of past contacts + redundancy
nLLA <- function( d, p, r ){
	#cat(r,p,"\n")
	pointLik <- by( d, d$cell, function(d){
		status <- tail(d$event.type,1)
		ncalc <- nrow(d)-1
		p.tot <- 0
		t <- d$time
		t <- head(tail( t, 1 ) - t,-1)

		if( status==1 ){
			p.tot <- sum( dgeom(0:(ncalc-1),p)*dexp(t,r) )
		} else {
			p.tot <- sum( dgeom(0:(ncalc-1),p)*(1-pexp(t,r)) )
			p.tot <- p.tot + 1 - pgeom(ncalc-1, p)
		}
		-log(p.tot)
	} ) 
	sum( pointLik )
}

fitA <- function(d){
	require(stats4)
	mle( function(p,r) nLLA(d,p,r), start=list( p=.2, r=.2 ), nobs=nrow(d),
		method = "L-BFGS-B", lower=list(p=.Machine$double.eps,r=.Machine$double.eps) )
}


#fa <- fitA()

reduce.d <- function( d ){ 
	require(stats4)
	fa <- fitA(d)
	p <- coef(fa)[['p']]
	r <- coef(fa)[['r']]
	Reduce( rbind, by( d, d$cell, function(d){
		status <- tail(d$event.type,1)
		ncalc <- nrow(d)-1
		t <- diff(tail(d$time,2))
		t <- d$time
		t <- head(tail( t, 1 ) - t,-1)
		if( status==1 ){
			d$V6 <- c(log(dgeom(0:(ncalc-1),p)*dexp( t, r )),-Inf)
		} else {
			d$V6 <- log(c(dgeom(0:(ncalc-1),p)*(1-pexp( t, r )),1-pgeom(ncalc-1,p)))
		}
		mi <- which.max( d$V6 )
		keep <- 1:nrow(d)
		keep <- !(d$event.type==2 & keep>mi)
		d$V7 <- as.integer( !keep )
		if( status == 1 ){
			d$V8 <- d$V6==max(d$V6)
		} else {
			d$V8 <- rep(0,nrow(d))
		}
		d
	} ) )
}

sim.h0 <- function( d, N=63 ){
	require(stats4)
	fa <- fitA(d)
	p <- coef(fa)[['p']]
	r <- coef(fa)[['r']]
	wt <- Reduce( c, by( d, d$cell, function(d) diff( d$time[d$event.type==2] ) ) )
	ct <- d[d$event.type<=0,"time"]
	#print( ct )
	o <- data.frame( cell=c(), cell_id=c(), time=c(), event.type=c() )
	for( i in seq_len( N ) ){
		id <- rep(i,length(t))
		ncalc <- 1+rgeom( 1,p )
		tim <- tail(diffinv(sample(wt,ncalc,replace=TRUE)),-1)
		last.t <- tail(tim,1)
		death.at <- rexp( r )+last.t
		censor.at <- sample( ct, 1 )
		if( death.at < censor.at ){
			# generate observed death
			t.redundant <- tail(diffinv(sample(wt,10*ncalc,replace=TRUE)),-1)+last.t
			t.redundant <- t.redundant[t.redundant < last.t]
			tim <- c( tim, t.redundant )
			o <- rbind( o, data.frame( cell=id, cell_id=id, time=c(tim,death.at),
				event.type=c(rep(2,length(tim)),1) ) )
		} else {
			# generate censored observation 
			tim <- tim[tim<=censor.at]
			o <- rbind( o, data.frame( cell=id, cell_id=id, time=c(tim,censor.at), 
				event.type=c(rep(2,length(tim)),0) ) )
		}
	}
	o
}

load.data <- function( reduce=FALSE, sim=FALSE ){
	d <- load.data.raw()
	if( reduce ){
		d <- reduce.d(d)
		write.table( d, "results/calcium-trace-reduced.txt",
			row.names=FALSE, col.names=FALSE )
		d <- d[d$V7==0,]
	}
	if( sim==TRUE ){
		d <- sim.h0( d )
	}
	time <- as.numeric(by( d$time, d$cell, function(d) diff(tail(d,2)) ))
	event <- as.numeric(by( d$event.type, d$cell, function(d) tail(d,1) ))
	previous.events <- as.numeric(by( d$event.type, d$cell, function(d) length(d)-2 ))
	# we don't distinguish between censoring by end of video and
	# censoring by emigration
	event[event==-1] <- 0

	data.frame( time=time, event=event, previous.events=previous.events )
}




s.to.df <- function( nps ){
		df.o <- data.frame( l=rep( names(nps$strata), nps$strata ) )
		df.o$t <- signif(nps$time)
		df.o$s <- signif(nps$surv)
		df.o
}

do.plot <- function( d, main, fname ){
	library( survival )
	#print( d )
	with( d ,{
		plot( survfit( Surv(time[previous.events==0], event[previous.events==0]) ~ 1 ),
			conf.int=FALSE, xlab="time after last Ca2+ hit (h)", xlim=c(0,24), 
			ylab="survival", main=main, mark.time=TRUE, col="black" )
		lines( survfit( Surv(time[previous.events==1], event[previous.events==1]) ~ 1 ),
			conf.int=FALSE, col="green", mark.time=TRUE )
		lines( survfit( Surv(time[previous.events==2], event[previous.events==2]) ~ 1 ),
			conf.int=FALSE, col="blue", mark.time=TRUE )
		lines( survfit( Surv(time[previous.events>2], event[previous.events>2]) ~ 1 ),
			conf.int=FALSE, col="red", mark.time=TRUE )
	} )
}

open.pdf <- function(fname="plots/figure-5b-km.R"){
	pdf(fname, width=4, height=4, useDingbats=FALSE)
	par( mar=c(4,4,1,1), bty="l", font.main=1, cex.main=1 )
}

open.pdf("plots/figure-5b-km.pdf")
do.plot( load.data(), "raw data", "results/raw" )
dev.off()

open.pdf("plots/figure-s4b-km-redundant-contacts-removed.pdf")
do.plot( load.data( TRUE ), "after ML removal of redundant contacts", "results/ml" )
dev.off()


open.pdf("plots/figure-5c-km-simulated.pdf")
do.plot( load.data( TRUE, sim=TRUE ), "simulated from H_0", "results/sim" )
dev.off()

