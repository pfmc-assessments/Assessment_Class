# FISH 576 & 577: Applied Stock Assessment I & II

This is a repository for the 2025 applied stock assessment course at the University of Washington.

## How this repository works:

This repository contains lectures and example code for both courses. [Something about it all being there or adding it as we go]

Please ensure you have a github handle and are subscribed to announcements under [discussions](https://github.com/okenk/Assessment_Class/discussions/categories/announcements). Github announcements will be our primary way of communicating with you as it allows students who drop the course to opt out, and students who start the course late to catch up on old communications. **Either check the github announcements regularly or ensure you are getting email alerts for new announcements.**

Feel free to also use discussions to communicate with your fellow classmates and the instructors!

## Instructors:

-   Vladlena Gertseva, NWFSC ([vladlena.gertseva\@noaa.gov](mailto:vladlena.gertseva@noaa.gov), [\@gertsevv](https://github.com/gertsevv))
-   Kiva Oken, NWFSC ([kiva.oken\@noaa.gov](mailto:kiva.oken@noaa.gov), [\@okenk](https://github.com/okenk))
-   Ian Taylor, NWFSC ([ian.taylor\@noaa.gov](mailto:ian.taylor@noaa.gov), [\@iantaylor-NOAA](https://github.com/iantaylor-NOAA))
-   Melissa Haltuch, AFSC ([melissa.haltuch\@noaa.gov](mailto:melissa.haltuch@noaa.gov)
-   Owen Hamel, NWFSC ([owen.hamel\@noaa.gov](mailto:owen.hamel@noaa.gov), [\@owenhamel](https://github.com/owenhamel))

## Overview:

Applied Stock Assessment I and II are a two-quarter applied stock assessment series offered in collaboration with stock assessment scientists in the Fishery Resource Analysis and Monitoring Division at the Northwest Fisheries Science Center. The first course provides a review of population dynamic modeling basics and stock assessment data types, and then focuses on the details of processing fishery and survey data for use in stock assessment and running Stock Synthesis (SS3) stock assessment models. The second course focuses on developing, running, evaluating and documenting the base and sensitivity Stock Synthesis (SS3) models required for a stock assessment update submitted to the Pacific Fishery Management Council (PFMC) for use in management.

## FISH 576: Applied Stock Assessment I

Students will work as a team to:

1.  Review assessment documents, STAR reports, and identify new literature

2.  Work up data for the update assessments

3.  Update data as each data source is finalized

4.  Begin work on producing an update stock assessment that involves updating and adding recent data from all data sources used in the previously reviewed stock assessment adopted for management.

### Learning goals:

Upon successful completion of the course, students will be able to:

1.  Evaluate and process length and age composition data and fishery-independent indices.

2.  Evaluate and process survey index data.

3.  Run existing Stock Synthesis models and replace or extend data in input files for catch, indices, composition, discard, and environmental data.

4.  Understand basic modeling assumptions and when they might be violated.

### Schedule:

| **Week** | **Topics**                                                                                                                             |
|------------|-----------------------------------------------------------|
| 1        | Class overview and an introduction                                                                                                     |
|          | Overview of stock assessment process, PFMC website                                                                                     |
|          | Introduction to git/github and class resources                                                                                         |
|          | Overview of data sources,non-disclosure forms, data preparation tasks for update assessments                                           |
|          | Review update assessment TORs                                                                                                          |
| 2        | WCGBTS/triennial/juvenile rockfish survey background and index standardization                                                         |
|          | Survey compositional data                                                                                                              |
|          | WCGBTS compositional data preparation, overview of nwfscSurvey code                                                                    |
| 3        | Fishery landings, PacFIN overview                                                                                                      |
|          | Fishery discards, WCGOP overview                                                                                                       |
| 4        | Fishery retained and discarded ages and lengths                                                                                        |
|          | PacFIN biological data processing and preparation using PacFIN.Utilities                                                               |
| 5        | Biological data and parameters estimated outside the model – Weight-Length, Maturity, Fecundity, Sex Ratios, Ageing precision and bias |
| 6        | Population modeling and assumptions, Integrated analysis,                                                                              |
|          | SS3 introduction, overview of SS3 input files                                                                                          |
|          | R4SS                                                                                                                                   |
| 7        | Modeling parameters in SS3                                                                                                             |
|          | Natural mortality and growth in SS3                                                                                                    |
|          | Initial conditions and fishing mortality in SS3                                                                                        |
|          | Recruitment, catchability and selectivity in SS3                                                                                       |
| 8        | Working with SS3 input files: formatting, processing and debugging                                                                     |
| 9        | Finalize data preparation                                                                                                              |
|          | Start running models                                                                                                                   |
| 10       | Data weighting in SS3                                                                                                                  |
|          | Formulate a proposed base model                                                                                                        |

## FISH 577: Applied Stock Assessment II

Students will work as a team to:

1.  Run the update stock assessment model and sensitivities.

2.  Build decision tables to provide catch advice to managers.

3.  Report results through a document and presentations.

Students will focus on running the stock assessment model under different configurations, ensuring model convergence, interpreting and comparing results and underlying assumptions, and documenting the stock assessment update in a document and presentations before the PFMC’s Scientific and Statistical Committee’s Groundfish Subcommittee. The PFMC’s Scientific and Statistical Committee and its Groundfish Subcommittee will review the update stock assessment.

### Learning goals:

Upon successful completion of the course, students will be able to:

1.  Read, understand and revise Stock Synthesis input files.

2.  Modify Stock Synthesis input files to reflect sensitivity runs, and run them.

3.  Run likelihood profiles in Stock Synthesis

4.  Produce forecasts using Stock Synthesis.

5.  Run R4SS software to produce plots from Stock Synthesis output.

6.  Write an updated stock assessment report.

7.  Present stock assessment results to the PFMC.

### Schedule:

| **Week** | **Topics**                                                                                                                                                                                                                                                                                                                                                                                                                     |
|------------|-----------------------------------------------------------|
| 1        | Review data topics from Applied Stock Assessment I                                                                                                                                                                                                                                                                                                                                                                             |
|          | Documentation for update stock assessments                                                                                                                                                                                                                                                                                                                                                                                     |
|          | Overview of R markdown document preparation                                                                                                                                                                                                                                                                                                                                                                                    |
| 2        | Model bridging: Run a set of models transitioning from the previous assessment model to the current by updating one new piece of data at a time and running the model. Plot comparisons of the previous assessment outputs with each subsequent model runs until each dataset has been updated. Produce a full set of r4ss output for the fully updated model.                                                                 |
| 3        | Risk neutrality and the science/management interface                                                                                                                                                                                                                                                                                                                                                                           |
|          | Retrospective analysis: Complete 5-year retrospective runs.                                                                                                                                                                                                                                                                                                                                                                    |
|          | Model convergence diagnostics, jittering.                                                                                                                                                                                                                                                                                                                                                                                      |
| 4        | Sensitivity analysis: Produce a set of model runs that include 1) model sensitivities from the last full (and any subsequent update) assessment, 2) any issues noted in the STAR or SSC reports, 3) runs that you are interested in completing, including any sensitivities that seem important given changes in parameter estimates in the update base model compared to the last model. Plot results against the base model. |
| 5        | Likelihood profiles: Complete likelihood profiles that were included in the last full assessment.                                                                                                                                                                                                                                                                                                                              |
| 6        | Management history/changes in management                                                                                                                                                                                                                                                                                                                                                                                       |
|          | Harvest Projections: Complete harvest projections and decision tables as included in the last full assessment.                                                                                                                                                                                                                                                                                                                 |
| 7        | Update stock assessment document preparation - NMFS internal review deadline (May19)                                                                                                                                                                                                                                                                                                                                           |
| 8        | Update stock assessment document preparation                                                                                                                                                                                                                                                                                                                                                                                   |
| 9        | Update stock assessment document preparation - PFMC briefing book deadline for documents (May 26)                                                                                                                                                                                                                                                                                                                              |
| 10       | Prepare and practice presentations for the PFMC SSC                                                                                                                                                                                                                                                                                                                                                                            |
|          | SSC GF review meeting (June 9-11)                                                                                                                                                                                                                                                                                                                                                                                              |
