Microperimetry and face recognition
======

These are R scripts to reproduce the analyses in Wallis, Taylor, Wallis, Jackson & Bex (under review).

These files are released under the GPL-3 License. I don't offer any support. You are welcome to send me an email but I make no claim that I will (be able to) help, particularly if you are attempting to use a different platform / setup than me (I'm running on OSX 10.8 with RStudio).

**IMPORTANT**: If you adapt any part of these scripts for your own academic work, please help keep me employed by citing the following paper:

Wallis, T.S.A., Taylor, C.P.T., Wallis, J., Jackson, M.L. and Bex, P.J. (under review). Characterisation of field loss based on microperimetry is predictive of face recognition difficulties.

Subdirectories in this repository
-----------------
The main file is manuscript.Rnw in the parent directory.
Other directories are:
  * /data/ contains all the raw data files from the experiment.
  * /funs/ contains all scripts and functions for the analysis.
  * /output/ contains files output by analysis. Notably, an R data file used for the analysis (output by /funs/import_data.R).
  * /figs/ contains all the figures in the manuscript.
  * /ideal_observer/ contains some legacy code used to produce the ideal observer analysis in the Supplementary Materials.


Data organisation
-----------------
Very raw data is contained in the /data/ subdirectory. If you're looking for a nice file to import all the data into your analysis program of choice, use /output/data.txt, which is a delimited text file including column headers. Each row is a trial in the face recognition experiment.

Reproducing the analysis
-----------------
The analyses in the paper were performed in R, using various libraries.
The main file is "manuscript.Rnw". This is an R Sweave file (to be compiled using
knitr).
This includes both the text of the manuscript and calls to the analysis functions
in the /funs/ subdirectory.
By looking through the manuscript.Rnw file you can determine which scripts in the
/funs/ directory to run in which order.
Alternatively, just compile manuscript.Rnw using knitr, for example, through Rstudio (www.rstudio.org).
You may need to uncomment some source() functions, for example, to run the MCMC sampler.

I will post more detailed instructions for reproducing the analysis soon.

License
-----------------
Copyright 2013, Thomas Wallis.

This is free, open source software released under the GPL-3 License. See LICENSE.txt for details.

Use this software at your own risk.
