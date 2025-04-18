# Changes to control files

## blocks
* end year of blocks may need to be adjusted, especially if it is equal to the final year of the model (could be changed to final year of new model)
* in some cases, new blocks will need to be added within the existing block patterns (e.g. explore adding additional block from 2011-2024 for retention parameters in widow model to allow Catch Shares era to differ from early period of fishery, as discussed in 13 March 2025 lecture)

## mortality and growth parameters
* prior on M may need to be updated (see Week 7 lecture: https://github.com/pfmc-assessments/Assessment_Class/blob/main/Class%20lectures/W07_M_Natural_Mortality.pdf)
* fixed weight-length parameters can be updated

## stock-recruit parameters
* last year of main recr_devs should be shifted as discussed in recruitment lecture
* the 5 recruitment bias adjustment values should be adjusted as discussed in recruitment lecture

## selectivity parameters
* if there are parameters hitting bounds, these may need to be adjusted to improve convergence
* if there is evidence that selectivity or retention has changed in recent years the blocks could be revised as noted above

## data weighting
* model needs to be retuned, with method used in full assessment, and new weights added to the control file


