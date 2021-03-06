# Location of foci of mitochondrial dysfunction along fibre section perimeters

[R code](LineScan2.R) and [three](Combined_manual_foci.txt) [data](Focix10.txt) [files](LineScans.txt) for muscle fibre section perimeter profiling and statistical analysis from [Vincent et al. (2018)](https://doi.org/10.1002/ana.25288).

The script [LineScan2.R](LineScan2.R) reads in data from a tab-delimted text file [Combinded_manual_foci.txt](Combined_manual_foci.txt) in the current working directory.  It uses the data to locate perinuclear regions around fibre perimeters and compares locations with manually specified foci (regions of mitochondrial dysfunction) along the same perimeters.  The script calculates observed overlap and overlap predicted by probability theory, assuming that locations of foci and perinuclei are independent.  As output, the script generates multi-page pdf reports demonstrating the locations of foci and perinuclei as well as their overlap.  The script also generates summary reports (histogram and scatterplot) comparing observed and predicted extent of overlap across dozens of fibre perimeters.  

We conclude from this analysis that foci co-localise with nuclei significantly more often than we might expect by chance.
