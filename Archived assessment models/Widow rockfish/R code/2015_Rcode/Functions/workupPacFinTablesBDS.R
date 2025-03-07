workupPacFinTablesBDS <- function(bds_fish,age_temp,sp_cluster,all_cluster) {
#############################################################################
#
# Get the output from SqlPlus queries of the Pacfin Hake data,
# combine and massage.
#
# Andi Stephens, 2010:  From John Wallace's code (his includes
#                       the queries, but this doesn't work well
#                       on all networks).
# Made as a function by Allan Hicks, November 2011
#
#  Input files:  bds.fish (from query)
#                bds.age  (from query)
#                bds.sp.cluster (from query)
#                bds.allsp.cluster (from query)
#
#  Output:       bds.fish.r.out.csv
#
#############################################################################


    #############################################################################
    #
    # Create unique record ids
    bds_fish$KEY <- paste(bds_fish$SAMPLE_YEAR, bds_fish$SOURCE_AGID, bds_fish$SAMPLE_NO, 
                          bds_fish$CLUSTER_NO, bds_fish$FISH_NO)
    age_temp$KEY <- paste(age_temp$SAMPLE_YEAR, age_temp$SOURCE_AGID, age_temp$SAMPLE_NO,
                      age_temp$CLUSTER_NO, age_temp$FISH_NO)
    sp_cluster$KEY <- paste(sp_cluster$SAMPLE_YEAR, sp_cluster$SOURCE_AGID, sp_cluster$SAMPLE_NO,
                      sp_cluster$CLUSTER_NO)
    all_cluster$KEY <- paste(all_cluster$SAMPLE_YEAR, all_cluster$SOURCE_AGID, all_cluster$SAMPLE_NO,
                      all_cluster$CLUSTER_NO)

    #############################################################################
    #
    # Check for ages, adjust bds_fish
    
    if(nrow(age_temp) == 0) {
      bds_fish$AGE_STRUCT_AGCODE <- ""
      bds_fish$AGE_METHOD <- ""
      bds_fish$AGE_READABILITY <- NA
      bds_fish$AGED_BY <- ""
      bds_fish$DATE_AGED <- as.POSIXct(NA)
      bds_fish$age1 <- NA
      bds_fish$age2 <- NA
      bds_fish$age3 <- NA    
    } else {
      # Sort age_temp including AGENUM to insure that AGED_BY and DATE_AGED are associated with the
      # first age (age1).
      col_names = c("SAMPLE_YEAR", "SOURCE_AGID", "SAMPLE_NO", "CLUSTER_NO", "FISH_NO", "AGENUM")
      age_temp <- age_temp[do.call(order, age_temp[ , col_names ]), ]

      # Remove duplicates, AGE_YEARS column
      bds_age <- age_temp[!duplicated(age_temp$KEY), ]
      bds_age$AGE_YEARS <- NULL

      # If there is a value for AGE_YEARS but not for AGENUM, replace AGE_YEARS with "-99"
      age_temp$AGE_YEARS[is.na(age_temp$AGE_YEARS) & !is.na(age_temp$AGENUM)] <- -99

      # Get AGE_YEARS for each record for each AGENUM, appending new columns
      # This lines up any multiple ages into columns so that one line = one fish
      for (i in unique(age_temp$AGENUM)) {
        bds_age <- cbind(bds_age, find.matching.rows(bds_age, age_temp[age_temp$AGENUM %in% i,],
                                        "KEY", "KEY", "AGE_YEARS"))
        colnames(bds_age)[ncol(bds_age)] <- paste("age", i, sep="")
      } # End for i
    
      # Clean up missing records
      if(is.null(bds_age$age2)) bds_age$age2 <- NA
      if(is.null(bds_age$age3)) bds_age$age3 <- NA
    
      bds_age$age1[is.na(bds_age$age1)] <- 0
      bds_age$age2[is.na(bds_age$age2)] <- 0
      bds_age$age3[is.na(bds_age$age3)] <- 0
    
      bds_age$age1[bds_age$age1 %in% -99] <- 0
      bds_age$age2[bds_age$age2 %in% -99] <- NA
      bds_age$age3[bds_age$age3 %in% -99] <- NA
       
      # Combine BDS_FISH and BDS_AGE
      col_names = colnames(bds_age)[grep("age", colnames(bds_age))]
      col_names = c("AGE_STRUCT_AGCODE", "AGE_METHOD", "AGE_READABILITY",
                    "AGED_BY", "DATE_AGED", col_names)
    
      bds_fish <- cbind(bds_fish, find.matching.rows(bds_fish, bds_age, "KEY", "KEY", col_names))
    
    } # End if-else
      
    #############################################################################
    # Now to clusters
    
    # Remove duplicates  
    sp_cluster <- sp_cluster[!duplicated(sp_cluster$KEY),]
    all_cluster <- all_cluster[!duplicated(all_cluster$KEY),]
    
    # The code below selects all clusters in a sample (regardless of species) and
    # then sums the cluster weight.  This is necesary only when there is a chance
    # of clusters that did not contain the target species.  The problem only seems
    # to occur in CA where the total weight of all clusters is not reported.           
    all_cluster.agg <- aggregate(list(all_cluster_sum = all_cluster$CLUSTER_WGT), 
                                   list(SAMPLE_NO = all_cluster$SAMPLE_NO), sum)
    sp_cluster <- cbind(sp_cluster, find.matching.rows(sp_cluster, all_cluster.agg,
                                            "SAMPLE_NO", "SAMPLE_NO", "all_cluster_sum"))
    
    #############################################################################
    # Combine SP_CLUSTER with BDS_FISH (which already has BDS_AGE if AGE = T)
    col_names = c("all_cluster_sum", "SPECIES_WGT", "CLUSTER_WGT", "FRAME_CLWT", "ADJ_CLWT")
                        
    # Reset bds_fish KEY to the same fields as sp_cluster KEY    
    bds_fish$KEY <- paste(bds_fish$SAMPLE_YEAR, bds_fish$SOURCE_AGID, bds_fish$SAMPLE_NO, 
                          bds_fish$CLUSTER_NO)
    bds_fish <- cbind(bds_fish, find.matching.rows(bds_fish, sp_cluster, "KEY", "KEY", col_names))
    
    #############################################################################
    # Duplicate all the records with frequency > 1 from Oregon
    bds_fish <- bds_fish[rep(1:nrow(bds_fish), bds_fish$FREQ),]
    
    #############################################################################
    # Cleanup
    bds_fish$KEY <- NULL
    
    # Return result
    return(bds_fish)
}
