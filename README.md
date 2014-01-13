Microperimetry and face recognition
======

These are R scripts to reproduce the analyses in Wallis, Taylor, Wallis, Jackson & Bex (2014).

These files are released under the GPL-3 License. I don't offer any support. You are welcome to send me an email but I make no claim that I will (be able to) help, particularly if you are attempting to use a different platform or setup than mine (I'm running on OSX 10.8 with RStudio).

**IMPORTANT**: If you adapt any part of these scripts for your own academic work, please help keep me employed by citing the following paper:

Wallis, T.S.A., Taylor, C.P.T., Wallis, J., Jackson, M.L. and Bex, P.J. (2014). Characterisation of field loss based on microperimetry is predictive of face recognition difficulties. _Investigative Ophthalmology and Visual Science, 55_(1): 142–153.

Subdirectories in this repository
-----------------
The main file is manuscript.Rnw in the parent directory.
Other directories are:
  * /data/ contains all the raw data files from the experiment.
  * /funs/ contains all scripts and functions for the analysis.
  * /output/ contains files output by analysis. Notably, an R data file used for the analysis (output by /funs/import\_data.R).
  * /figs/ contains all the figures in the manuscript.
  * /ideal\_observer/ contains some legacy code used to produce the ideal observer analysis in the Supplementary Materials.


Data organisation
-----------------
Very raw data is contained in the /data/ subdirectory. If you're looking for a nice file to import all the data into your analysis program of choice, use /output/data.txt, which is a delimited text file including column headers. Each row is a trial in the face recognition experiment.

Reproducing the analysis
-----------------
The analyses in the paper were performed in R, using various libraries.
The main file is "manuscript\_full.Rnw". This is an R Sweave file (to be compiled using
knitr).
The manuscript\_full.Rnw file will produce the full manuscript, including tables and the supplementary file.
This includes both the text of the manuscript and calls to the analysis functions
in the /funs/ subdirectory.
By looking through the manuscript.Rnw file you can determine which scripts in the
/funs/ directory to run in which order.
Alternatively, just compile manuscript.Rnw using knitr, for example, through Rstudio (www.rstudio.org).
You may need to uncomment some source() functions, for example, to run the MCMC sampler.

Dependencies
-----------------
R packages that you will need to install to reproduce the analyses in full are
  * Rstan (available from www.mc-stan.org)
  * plyr (on CRAN)
  * ggplot2 (on CRAN)
  * reshape2 (on CRAN)
  * wutils (available from my github [here](https://github.com/tomwallis/wutils))
  * psybayes (available from my github [here](https://github.com/tomwallis/psybayes))


Stan versions for the sampling in the paper
-----------------

The MCMC samples for the main analysis (Figures 3 -- 5), and for the model fit to only block 1, were generated using Stan version 1.3. The samples for the correlations (Table 3) were generated using Stan version 2.0.0.

The paper was accepted on November 4, 2013. On 27 December 2013 Stan 2.1.0 was released, and included a major bug fix for Stan versions 2.0.0 and 2.0.1:

> Major Bug in 2.0.0, 2.0.1
> ------------------------------
> Stan 2.0.0 and Stan 2.0.1 introduced a bug in the implementation
> of the NUTS criterion that led to poor tail exploration and
> thus biased the posterior uncertainty downward.  There was no
> bug in NUTS in Stan 1.3 or earlier, and 2.1 has been extensively tested
> and tests put in place so this problem will not recur.

The main analyses are not affected by this bug because they were run under Stan version 1.3. I have re-run the correlations (Table 3) with Stan 2.1.0 and the substantive results do not change. Fixation stability and visual acuity are both correlated with percent correct; no other correlations were credibly different to zero.

License
-----------------
Copyright 2013, Thomas Wallis.

This is free, open source software released under the GPL-3 License. See LICENSE.txt for details.

Use this software at your own risk.
