                             
READ ME (please)

--------------------------------------------------------------------
OVERVIEW

Data and R code to reproduce the analyses in the manuscript:

“Novel parasite invasion destabilizes species coexistence” 

Which is a paper about the impact of an invasive parasitic nematode on the population densities and demographic rates of killifish and guppies living in two community types in streams in a beautiful rainforest in Trinidad.

By Tomos Potter*, Ryan S. Mohammed, Joshua Goldberg, David N Reznick, Joseph Travis, and Ronald D Bassar

* you can contact me here: tomos.potter@protonmail.com 

IMPORTANT: To run this code, you will first need to download and install Program Mark, which is freely available here: http://www.phidot.org/software/mark/downloads/

You won't have to open Program Mark yourself - all work is done via R, using the RMark() package to interface with the external software.

Note that the scripts use relative pathways to access data, call other scripts, save output etc, so please keep the same file names and structure as in this repository if you want everything to run nice and smooth.

--------------------------------------------------------------------
CONTENTS

The main repository contains the following folders and files (descriptions in parentheses):

Files:

-killifish-guppy-parasite-mark-recapture.Rproj
(The RStudio project file for this project. Click on it to launch the project in RStudio)

- READ-ME.txt 
(yes, this very thing that you are reading right now)

- MAIN-SCRIPT.R 
(This script runs the full analysis by calling all the other scripts, loading the data, running the models, saving the results to file, plotting the figures etc. ***This is the only script you need to run to generate the results***)

Folders and files within:

- DATA 
(contains the 4 data files used in the analyses, plus 2 images used in the figures):

 1. killifish-mr-data.csv (mark-recapture data for killifish)
 2. guppy-mr-data.csv (mark recapture data for guppies)
 3. dissection-data.csv (data from dissected fish)
 4. stream-data.csv (data on stream area to calculate population densities)
 5. killifish_silhouette.png (killifish image used in figures)
 6. guppy_silhouette.png (guppy image used in figures)

- OUTPUT 
(this is where model outputs, tables of results, and figures are saved to. I have included all of the output here so you can examine them if you like without having to run the models (i.e. if you don't have the time or don't want to download Program Mark). There is a lot of stream-specific output, I.e. either from Caigual (CA) or Taylor (TY). If you run the main script, it will save over the files that are present here. 

There are two sub-folders: figures, and results.
	
	- figures

(This is where the figures in the manuscript are saved to)

 1. Figure-1.png
 2. Figure-2.png
 3. Figure-3.png
 4. Figure-4.png
 5. Figure-5.png
 


	- results

(here are 20 (!) output objects: 6 model objects (with extension .Rda), and 14 .csv files of results / parameter estimates)

 1. AICCA.csv			(killifish POPAN AIC table for CA)
 2. AICCAsize.csv 		(killifish Multistrata AIC table for CA)
 3. AICTY.csv			(killifish POPAN AIC table for TY)
 4. AICTYsize.csv 		(killifish Multistrata AIC table for TY)
 5. G-AICCA.csv			(guppy POPAN AIC table for CA)
 6. G-AICTY.csv			(guppy POPAN AIC table for TY)
 7. G-POPAN-parameters.csv	(results: guppy POPAN parameters)
 8. G-POPAN-model-CA.Rda	(all guppy POPAN models for CA)
 9. G-POPAN-model-TY.Rda	(all guppy POPAN models for TY)
10. guppy-infection-rate-dissection.csv (results from guppy dissection)
11. guppy-infection-rate-focal.csv (guppy infection rate in CA and TY)
12. killifish-infection-rate-dissection.csv (results from killifish dissection)
13. killifish-infection-rate-focal.csv (killifish infection rate in CA and TY)
14. MPMs-toplot.csv		(results from matrix projection models)
15. Multistrata-models-CA.Rda	(all killifish Multistrata models for CA)
16. Multistrata-models-TY.Rda	(all killifish Multistrata models for TY)
17. Multistrata-parameters.csv	(results: killifish Multistrata parameters)
18. POPAN-models-CA.Rda		(all killifish POPAN models for CA)
19. POPAN-models-TY-.Rda	(all killifish POPAN models for TY)
20. POPAN-parameters.csv	(results: killifish POPAN parameters


- SUB-SCRIPTS 
(contains 6 R scripts, which are called by MAIN-SCRIPT.R)

 1. POPAN-killifish.R (runs the POPAN models for killifish)
 2. POPAN-guppies.R (runs the POPAN models for guppies)
 3. Multistrata-killifish.R (runs the Multistrata models for killifish)
 4. MPMs.R (runs the matrix projection models for killifish)
 5. Parasite-prevalence.R (runs analysis of parasite prevalence in both species)
 6. Figures.R (unsurprisingly, this one makes the figures)

--------------------------------------------------------------------
END
