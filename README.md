# Calcium Event Analysis Code

This repository contains the data and code for the analysis of calcium signalling events performed for the pape "T cell-mediated additive cytotoxicity to effectively eliminate cancer cells" by Weigelin et al., to be published in Nature Communications (a preprint is available here: https://www.biorxiv.org/content/10.1101/2020.04.22.054718v1). 

The person responsible for these analyses was Johannes Textor (current contact information is available at http://johannes-textor.name). If you have questions or comments on these analyses, please contact Johannes Textor; otherwise, contact the lead and corresponding author, Bettina Weigelin (contact information is available through the link above). 

Specifically, this repository includes code for the following figure panels:

 * Figure 5b,c(TODO),d
 * Supplementary Figure 4a,b(TODO)

If you are looking for data or code underlying the other figures and panels, please contact Bettina Weigelin. 

There are also two scripts, prefixed with "extra", that run analyses that were requested by the reviewers and included in the rebuttal letter, but not in the final paper. We provide these for transparency only. 

# Running the code

To run everything, you should just need to type `make`. The generated figures will be stored in the folder `plots/`. These figures are named correspondingly to how they appear in the paper. If you want to re-do everything, run `make -B` or run `make clean ; make`.

Those who don't have Python or make can also simply run the individual R scripts that make each figure panel. The only dependency for these scripts is the file `data/calcium-trace.txt`, which is also included in this repository. 

# Data

The original data, as recorded during the experiment, is available as an Excel sheet in the folder "data/". For the purpose of this analysis, the Python script `import-excel-data.py` converts the data into a more convenient form, a space-separated text file `data/calcium-trace.txt` containing four columns:

"cell" -- a consecutive identifier given to each cell
"cell_id" -- the original ID of the cell as used in the Excel file (not used in the analysis)
"time" -- the time at which the event took place
"event.type" -- the type of event:

-1 -- cell leaves imaging region 
0  -- end of movie
1  -- cell dies
2  -- Ca2+ flux

"ctl" -- for Ca2+ events, the number of the associated CTL that this event came from (as judged by the person who analyzed the movies). 

