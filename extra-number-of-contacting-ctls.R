

## Extra analysis used in rebuttal letter only, but not in the paper:

## Show how many CTL contacts happened before each apoptosis.

d <- read.table("data/calcium-trace.txt", header=TRUE)[,-2]
d$event.type[d$event.type==-1] <- 0

# Extract only those cells that have had at least 1 Ca2+ event
#nr <- c( by(d, d$cell, nrow) )
#d <- d[d$cell %in% names(which(nr>1)),]

# cell fate: 1 = dead, 0 = alive
fate <- c( by(d$event.type, d$cell, function(x) tail(x,1)) )

# focus on dead cells only
d <- d[d$cell %in% names(fate[fate==1]),]

n.ctl <- function(tlag){
	tabulate( 1+by( d, d$cell, function(x){
		x <- x[x$time >= tail(x$time,1)-tlag,]
		length( unique( x[x$event.type==2,]$ctl ) )
	} ), 6 )
}

x <- sapply( c(seq(0,6,by=0.5),Inf), n.ctl )

colnames( x ) <- c(paste0( seq(0,360,by=30), " min"), "overall") 

rownames( x ) <- paste0( seq_len(nrow(x))-1, " CTL" )

x <- t(x)

writexl::write_xlsx( cbind( `time before death`=rownames(x), 
	as.data.frame(x)), "data/contact-profile.xlsx" )
