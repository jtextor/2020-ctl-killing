

d <- read.table("data/calcium-trace.txt", header=TRUE)[,-2]
d$event.type[d$event.type==-1] <- 0

plt <- function(d){
	n.ca.hits <- c(by( d, d$cell, function(x)nrow(x)-1 ))
	cells <- unique(d$cell)

	final.events.only <- which( d$event.type <= 1 )

	end.times <- d$time[final.events.only]
	fate <- d$event.type[final.events.only]
	w.ca.hits <- which(d$event.type==2)
	col.ctl <- unlist( by( d$ctl, d$cell, function(x) c(1,1+cumsum(diff(x))) ) )
	ns <- seq_along(end.times)

	plot(NA, xlim=c(0,25), xaxt="n", ylim=c(1,length(end.times)),
			yaxt="n", bty="n", xlab="time (h)")
	axis(1,at=seq(0,24,by=6))
	segments(0,ns,end.times,ns)
	points(end.times,ns,pch=21-fate, cex=c(1.4,1)[2-fate])


	wcell <- as.integer(factor(d$cell))

	points(  d$time[w.ca.hits], wcell[w.ca.hits], pch=21, cex=.6, 
			bg=1+col.ctl[w.ca.hits] )
}

pdf("plots/figure-s4a-swimmer.pdf", width=9, height=9 )
par(mar=c(4,0.2,0.2,0.2),font.main=1, cex.main=1)
plt(d)
dev.off()

