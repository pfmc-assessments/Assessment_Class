# FISH 576 & 577: Applied Stock Assessment I & II

This is a repository for the 2027 applied stock assessment course at the University of Washington. 

## How this repository works:

This repository contains lectures, example code, and background reading materials for both courses. Materials will be added as the course progresses.

Please ensure you have a github handle and are subscribed to announcements under [discussions](https://github.com/okenk/Assessment_Class/discussions/categories/announcements). Github announcements will be our primary way of communicating with you as it allows students who drop the course to opt out, and students who start the course late to catch up on old communications. **Either check the github announcements regularly or ensure you are getting email alerts for new announcements.**

Feel free to also use discussions to communicate with your fellow classmates and the instructors!

## Instructors:

-   Dr. Vladlena (Vlada) Gertseva, NWFSC ([vladlena.gertseva\@noaa.gov](mailto:vladlena.gertseva@noaa.gov), [\@gertsevv](https://github.com/gertsevv))
-   Dr. Ian Taylor, NWFSC ([ian.taylor\@noaa.gov](mailto:ian.taylor@noaa.gov), [\@iantaylor-NOAA](https://github.com/iantaylor-NOAA))

## Overview:

Applied Stock Assessment I and II are a two-quarter applied stock assessment series offered in collaboration with stock assessment scientists in the Fishery Resource Analysis and Monitoring Division at the Northwest Fisheries Science Center. The first course provides a brief review of population dynamic modeling basics and stock assessment data types, and then focuses on the details of processing fishery and survey data to develop a stock assessment update using Stock Synthesis (SS3) modeling platform. The second course focuses on developing, running, evaluating and documenting the Stock Synthesis (SS3) models required for a stock assessment update submitted to the Pacific Fishery Management Council (PFMC) for use in management. This year, the course also introduces the Fisheries Integrated Modeling System (FIMS), a next-generation stock assessment framework, to run the FIMS model in parallel with SS3, to help with future transition to FIMS.  

## Course format:

This is a hybrid course, with in-person and virtual participation available. In person, we will meet at UW in FSH 105. For online login information, contact the instructors.

## Scheduled course times (subject to change based on participant schedules): 

Tuesdays 11:00am-11:50am and Thursdays 10:30am-11:20am

## FISH 576: Applied Stock Assessment I

Students will work as a team to:

1. Learn about the fisheries management system in the U.S. West Coast; review assessment documents, stock assessment review (STAR) reports, and identify new relevant literature.

2.  Become acquainted with Stock Synthesis and FIMS.

3.  Work up data for the update stock assessments for each data source.

4.  Begin producing an update assessment that involves updating and adding recent data from all data sources used in the previously reviewed stock assessment adopted for management.

### Learning goals:

Upon successful completion of the course, students will be able to:

1.  Evaluate and process fishery catch and survey index data..

2.  Evaluate and process length and age composition data..

3.  Run existing Stock Synthesis model and replace or extend data in input files for catch, indices, composition, discard, and environmental data.

4.  Explore current features and start running the FIMS model

5.  Understand basic modeling assumptions and when they might be violated.

### Schedule:

| **Week** | **Topics**                                                                                                                             |
|------------|-----------------------------------------------------------|
| 1        | Class overview and an introduction                                                                                                     |
|          | Overview of stock assessment process, PFMC website                                                                                     |
|          | Introduction to git/github and class resources                                                                                         |
|          | Overview of data sources,non-disclosure forms, data preparation tasks for update assessments                                           |
|          | Review update assessment TORs                                                                                                          |
| 2        | Population modeling and assumptions, Integrated analysis,                                                                              |
|          | SS3 introduction, overview of SS3 input files                                                                                          |
|          | Working with SS3 input files: formatting, processing and debugging                                                                     |
|          | R4SS                                                                                                                                   |
| 3        | Introduction to FIMS                                                                                                                   |
|          | Overview of FIMS current and future features, FIMS resources                                                                           |
| 4        | WCGBTS/triennial survey background and index standardization                                                                           |
|          | Survey compositional data                                                                                                              |
|          | WCGBTS compositional data preparation, overview of nwfscSurvey code                                                                    |
| 5        | Fishery landings, PacFIN overview                                                                                                      |
|          | Fishery discards, WCGOP overview                                                                                                       |
| 6        | Fishery retained and discarded ages and lengths                                                                                        |
|          | PacFIN biological data processing and preparation using pacfintools                                                                    |
| 7        | Biological data and parameters estimated outside the model – Weight-Length, Maturity, Fecundity, Sex Ratios, Ageing precision and bias |
|          | Recruitment index based on oceanographic data                                                                                          |
| 8        | Modeling parameters                                                                                                                    |
|          | Natural mortality and growth                                                                                                           |
|          | Initial conditions and fishing mortality                                                                                               |
|          | Recruitment, catchability and selectivity                                                                                              |
| 9        | Finalize data preparation                                                                                                              |
|          | Data weighting                                                                                                                         |
| 10       |  Formulate a proposed base model                                                                                                       |
| Final week    | Present proposed base model in SS3, and FIMS parallel model                                                                       |

## FISH 577: Applied Stock Assessment II

Students will work as a team to:

1.  Run the update stock assessment model and complete a set of required model diagnostics.

2.  Build decision tables to provide catch advice to managers.

3.  Communicate results through a written report and presentations.

Students will focus on running the stock assessment model under different configurations and underlying assumptions, ensuring model convergence, interpreting and comparing results, and documenting the stock assessment update in an assessment report and presentations before the PFMC’s Scientific and Statistical Committee’s Groundfish Subcommittee. The PFMC’s Scientific and Statistical Committee and its Groundfish subcommittee review the update stock assessment.

### Learning goals:

Upon successful completion of the course, students will be able to:

1.  Read, understand and modify Stock Synthesis and FIMS input files.

2.  Produce model results and plots from model outputs.

3.  Run model diagnostics, such as likelihood profiles, sensitivity and retrospective analyses.

4.  Develop assessment model forecasts, for use in management.

5.  Write a detailed stock assessment report.

6.  Present stock assessment results to the PFMC and potentially other stake holders.

### Schedule:

| **Week** | **Topics**                                                                                                                                                                                                                                                                                                                                                                                                                     |
|------------|-----------------------------------------------------------|
| 1        | Review data topics from Applied Stock Assessment I                                                                                                                            |
|          | Documentation for update stock assessments                                                                                                                                    |
|          | Overview of  asar document preparation process                                                                                                                                |
|          | Overview of diagnostics to be complete                                                                                                                                        |
| 2        | Model bridging: Run a set of models transitioning from the previous assessment model to the current by updating one new piece of data at a time and running the model. Plot comparisons of the previous assessment outputs with each subsequent model runs until each dataset has been updated. Produce a full set of R4SS output for the fully updated model.                                                                  |
| 3        | Model convergence diagnostics, jittering                                                                                                                                      |
|          | Retrospective analysis: Complete 5-year retrospective runs                                                                                                                    |
|          | Generate assessment report template and start writing the stock assessment document                                                                                           |
| 4        | Sensitivity analysis: Produce a set of model runs that include 1) model sensitivities from the last full (and any subsequent update) assessment, 2) any issues noted in the STAR or SSC reports, 3) runs that you are interested in completing, including any sensitivities that seem important given changes in parameter estimates in the update base model compared to the last model. Plot results against the base model. |
|          | Complete Introduction, Data and Model Description sections of the assessment document                                                                                         |
| 5        | Likelihood profiles: Complete likelihood profiles that were included in the last full assessment                                                                              |
|          | Complete Model Diagnostics section of the assessment document                                                                                                                 |
| 6        | Management history/changes in management                                                                                                                                      |
|          | Risk neutrality and the science/management interface                                                                                                                          |
|          | Harvest Projections: Complete harvest projections and decision tables as included in the last full assessment.                                                                |
|          | Complete Executive Summary of the assessment document                                                                                                                         |
| 7        | Complete stock assessment document for instructors’ review                                                                                                                    |
| 8        | Complete stock assessment document for NMFS/PFMC internal review deadline                                                                                                     |
| 9        | Complete stock assessment document for  PFMC briefing book deadline for the SSC review                                                                                        |
| 10       | Prepare and practice presentation for the PFMC’s SSC Groundfish Subcommittee review meeting                                                                                   |
| Final week   | Present update stock assessment to the PFMC’s Scientific and Statistical Committee’s Groundfish Subcommittee                                                              |
