.SECONDARY:
.DELETE_ON_ERROR:

all : data/calcium-trace.txt plots/figure-s4a-swimmer.pdf \
	plots/figure-5b-km.pdf plots/figure-5d-estimate-repair-time.pdf

data/calcium-trace.txt : import-excel-data.py data/ca2-data.xlsx
	python3 $^ >$@

plots/%.pdf : %.R data/calcium-trace.txt
	Rscript $<

clean:
	rm -f plots/* data/calcium-trace.txt
