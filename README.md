# [QFA](http://research.ncl.ac.uk/qfa/) screens to search for genetic interaction with a *bir1-17* mutation

Quantitative Fitness Analysis ([QFA](http://research.ncl.ac.uk/qfa/)) is a set of experimental and computational techniques for analysing the growth of microbial colonies arrayed on the surface of solid agar plates.  QFA is usually used to search for genetic interactions with a query mutation.  In this case, we look for genetic interactions with *bir1-17*, a temperature sensitive allele of *[BIR1](http://www.yeastgenome.org/locus/BIR1/overview)*.  

To find genetic interactions with *bir1-17* we carried out a screen of the yeast knockout collection crossed with a *bir1-17* query strain to give ~5,000 double mutant strains.  We also generated a second set of ~5,000 control strains by crossing the same knockout collection with a wild-type strain (labelled cSGA).

This repository contains four directories whose contents are described below.  An R script file (updateFiles.R) and this file (README.md) can be found in the root directory.

## BIR1-17

### ANALYSISOUT
Fitness estimates generated from raw growth curve data in IMAGELOGS directory by fitting the logistic model to data and by smoothing raw data, using the [QFA R package](http://qfa.r-forge.r-project.org/).

### IMAGELOGS
Raw [Colonyzer](http://research.ncl.ac.uk/colonyzer) image analysis output for four replicate *bir1-17* screens.

## cSGA

### ANALYSISOUT
Fitness estimates generated from raw growth curve data in IMAGELOGS directory by fitting the logistic model to data and by smoothing raw data, using the [QFA R package](http://qfa.r-forge.r-project.org/).

### IMAGELOGS
Raw [Colonyzer](http://research.ncl.ac.uk/colonyzer) image analysis output for four replicate cSGA screens.

## qfaDALBIRHISTORICAL
Files for building the [QFA visualisation tool](http://qfa.r-forge.r-project.org/visTool) for this project.

## GIS_output
Output report tables reporting Genetic Interaction Strengths (GIS.txt), fitness plots for visualising the evidence for genetic interaction (.pdf)

