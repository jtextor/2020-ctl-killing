.SECONDARY:
.DELETE_ON_ERROR:

all : data/calcium-trace.txt plots/swimmer.pdf plots/km.pdf

data/calcium-trace.txt : import-excel-data.py data/ca2-data.xlsx
	python3 $^ >$@

plots/%.pdf : %.R data/calcium-trace.txt
	Rscript $<
