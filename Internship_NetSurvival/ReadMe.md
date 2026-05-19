# Getting started with the project

This material is intended to help you become familiar with the concepts and tools needed for the internship.

You are not expected to understand everything in detail at this stage. The goal is to progressively build intuition.

## 1. Video lectures (start here)

I recommend starting with the following two video lectures, which provide a clear and intuitive introduction to the key concepts:

- "Non-parametric estimators of net survival". This is a lecture given by Maja Pohar-Perme introducing net survival and population-based indicators. 
- "Flexible parametric models for excess hazard and introduction to penalized framework. Part 1". Focus on the first talk given by Emmanuelle Dantony about excess mortality regression models.

You can access the material here:

https://sesstim.univ-amu.fr/fr/page/pedagogical-material-corsican-summer-school-2024


Password: corsican2024

This will give you access to a table with slides, videos, and practical sessions.

At this stage, I recommend focusing only on these two videos.
The platform contains many additional materials, but going through everything is not necessary and may be overwhelming at this point. You can come back to the other materials later if needed.

## 2. Core reading 

Please start with the following papers:

- PoharPermeEsteveRachet2016_NetSurvivalControversies.pdf: overview of net survival and related concepts.

- PoharPermePavlic2018_relsurvPackage.pdf: introduction to non-parametric estimation and R implementation.

- PoharPermeHendersonStare2009_RelativeSurvivalRegression.pdf: regression framework and methodological aspects. You also have a pdf containing the supplementary material.

These will give you the conceptual foundations of net survival.


## 3. Data generation (important for the project)

Understanding how data are generated is important for the simulation study.

- The extract from my PhD (in French) describes the data generation process, in particular how event times are simulated. This is the file PhDExtract_DataGeneration_ExcessHazard_FR.pdf.

Note that the video lectures provide a conceptual understanding of the excess mortality regression model, which underlies this generation process.

The key idea is that, at the individual level (i.e., conditional on covariates), the total hazard can be expressed as the sum of the expected hazard and the excess hazard.

## 4. Code

You will also find an initial R code base.

At this stage, the goal is simply to:

- run the code
- understand its structure
- identify the main components

You will later extend and adapt it when designing the simulation study.

## 5. First steps

Suggested first steps:

1. Watch the two recommended videos
2. Read the core papers
3. Run the code
   
Take notes and list questions: we will discuss your understanding and next steps together.