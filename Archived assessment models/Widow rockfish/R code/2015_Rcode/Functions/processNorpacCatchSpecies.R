processNorpacCatch <- function(ncatch,species=206,outfname=NULL) {
    ##########################
    #Process NORPAC domestic database extraction
    #Code developed by Andi Stephens in 2010
    #Reformatted by Allan Hicks in 2011
    #
    #  Domestic Norpac catch consists of sampled hauls and unsampled hauls.
    #  For sampled hauls, EXTRAPOLATED_WEIGHT is the catch.  For unsampled
    #  hauls, we get a monthly bycatch rate and apply it to reduce the
    #  OFFICIAL_TOTAL_CATCH.
    #
    #  Assumptions:
    #
    #    1.  EXTRAPOLATED_WEIGHT is in units 1,000 times those of the
    #        OFFICIAL_TOTAL_CATCH
    #    2.  All sampled hauls have an entry for SPECIES == 206, i.e., Hake.
    #    3.  Each unsampled haul is represented by one unique record.
    ##############################
    
    colnames(ncatch)[1] <- "HAULJOIN"
    ncatch$SPECIFICHAUL <- paste(format(ncatch$HAULJOIN,digits=19),ncatch$HAUL,sep="_")
    ncatch$Month <- as.numeric(substring(ncatch$HAUL_DATE,6,7))
    ncatch$Year <- as.numeric(substring(ncatch$HAUL_DATE,1,4))

    years <- sort(unique(ncatch$Year))

    allOut <- NULL 
    for (i in years) {
        x <- ncatch[ncatch$Year == i,]  #only the specific year
        x <- x[(!is.na(x$OFFICIAL_TOTAL_CATCH)) & x$OFFICIAL_TOTAL_CATCH != 0,]  #only observations with an official catch
        
        Nrec <- nrow(x)
        Nhauls <- length(unique(x$SPECIFICHAUL))
        
        #Adjust units so that the OTC and EXTRAPOLATED_WEIGHT are in same units.
        x$OFFICIAL_TOTAL_CATCH <- x$OFFICIAL_TOTAL_CATCH * 1000
        
        #Dummy for now
        x$ActualCatch <- x$EXTRAPOLATED_WEIGHT
        
        #Get a raw bycatch.  This will be aggregated per month.
        #x$ByCatch <- 0
        
        #Separate the unsampled hauls.  Only one record per haul if unsampled
        unsampled <- x[is.na(x$HAUL_SAMPLED_BY) | x$HAUL_SAMPLED_BY == 0,]
        Un.hauls <- length(unique(unsampled$SPECIFICHAUL))
        
        cat("Year", i, "\tHauls", Nhauls, "\tUnsampled", Un.hauls, "\n"); flush.console();
        
        #Can only calculate a rate factor from the HAKE catch. Bycatch is the official total weight minus the hake catch (sp=206)
        #Use only sampled hauls and hake to calculate bycatch
        x <- x[(!is.na(x$HAUL_SAMPLED_BY)) & x$HAUL_SAMPLED_BY != 0 & x$SPECIES == species,]
        x$ByCatch <- x$OFFICIAL_TOTAL_CATCH - x$EXTRAPOLATED_WEIGHT 
        
        if (nrow(unsampled) > 0) { 
            # Use the sampled hauls to generate a monthly bycatch rate.
        
            B.byMonth <- aggregate(x$ByCatch, list(x$Month), sum, na.rm=T)
            T.byMonth <- aggregate(x$OFFICIAL_TOTAL_CATCH, list(x$Month), sum, na.rm=T)
#print(dim(T.byMonth))
            names(B.byMonth) <- c("Month", "ByCatch")
            names(T.byMonth) <- c("Month", "Catch")
            months <- unlist(T.byMonth[1])
            names(months) <- NULL
        
            bycatch.rate <- B.byMonth[,2] / T.byMonth[,2]
        
            Fraction.of.total <- rep(NA,12) #for each month
            Fraction.of.total[months] <- 1-bycatch.rate
#print(Fraction.of.total)
            #Apply bycatch rate to unsampled hauls with total catch
            unsampled$ActualCatch <- unsampled$OFFICIAL_TOTAL_CATCH * Fraction.of.total[unsampled$Month] 
            unsampled$ByCatch <- unsampled$OFFICIAL_TOTAL_CATCH - unsampled$OFFICIAL_TOTAL_CATCH * Fraction.of.total[unsampled$Month] 
                
            # Put the dataset back together
            x <- rbind(x, unsampled)
        } # End if there are unsampled hauls.
        
        # Not output, but might want it...
        #Raw.Total <- aggregate(x$OFFICIAL_TOTAL_CATCH, list(x$Month), sum, na.rm=T)
        
        # Generate Monthly Catch and return it
        M.C <- aggregate(x$ActualCatch, list(x$Month), sum, na.rm=T)
        colnames(M.C) <- c("Month", "Catch")
        #M.C$Year <- i
        #M.C$Sector <- "Domestic"
        
        outdata <- data.frame(Sector="DomesticAtSea", Month= M.C$Month, Year=i, Catch.MT=M.C$Catch/1000)
        #colnames(outdata) = c("Month", "Year", "Catch", "Sector")
        allOut <- rbind(allOut,outdata)
        
        if(!is.null(outfname)) {
            if (i == years[1]) {
                #    write.table(file=outfname, outdata, col.names=T, row.names=F, sep=",", append = T)
                write.csv(file=outfname, outdata, col.names=T, row.names=F, append = F)
            } else {
                write.csv(file=outfname, outdata, col.names=F, row.names=F, sep=",", append = T)
            } # End if-else
        }
    } # End for i
    
    invisible(allOut)
}



processNorpacCatchDetail <- function(ncatch,outfname=NULL) {
    ##########################
    #Process NORPAC domestic database extraction
    #Code developed by Andi Stephens in 2010
    #Reformatted by Allan Hicks in 2011
    #
    #  Domestic Norpac catch consists of sampled hauls and unsampled hauls.
    #  For sampled hauls, EXTRAPOLATED_WEIGHT is the catch.  For unsampled
    #  hauls, we get a monthly bycatch rate and apply it to reduce the
    #  OFFICIAL_TOTAL_CATCH.
    #
    #  Assumptions:
    #
    #    1.  EXTRAPOLATED_WEIGHT is in units 1,000 times those of the
    #        OFFICIAL_TOTAL_CATCH
    #    2.  All sampled hauls have an entry for SPECIES == 206, i.e., Hake.
    #    3.  Each unsampled haul is represented by one unique record.
    ##############################
    
    colnames(ncatch)[1] <- "HAULJOIN"
    ncatch$SPECIFICHAUL <- paste(format(ncatch$HAULJOIN,digits=19),ncatch$HAUL,sep="_")
    ncatch$Day <- as.numeric(substring(ncatch$HAUL_DATE,9,10))
    ncatch$Month <- as.numeric(substring(ncatch$HAUL_DATE,6,7))
    ncatch$Year <- as.numeric(substring(ncatch$HAUL_DATE,1,4))

    years <- sort(unique(ncatch$Year))

    allOut <- NULL 
    for (i in years) {
        x <- ncatch[ncatch$Year == i,]  #only the specific year
        x <- x[(!is.na(x$OFFICIAL_TOTAL_CATCH)) & x$OFFICIAL_TOTAL_CATCH != 0,]  #only observations with an official catch
        
        Nrec <- nrow(x)
        Nhauls <- length(unique(x$SPECIFICHAUL))
        
        #Adjust units so that the OTC and EXTRAPOLATED_WEIGHT are in same units.
        x$OFFICIAL_TOTAL_CATCH <- x$OFFICIAL_TOTAL_CATCH * 1000
        
        #Dummy for now
        x$ActualCatch <- x$EXTRAPOLATED_WEIGHT
        
        #Get a raw bycatch.  This will be aggregated per month.
        #x$ByCatch <- 0
        
        #Separate the unsampled hauls.  Only one record per haul if unsampled
        unsampled <- x[is.na(x$HAUL_SAMPLED_BY) | x$HAUL_SAMPLED_BY == 0,]
        Un.hauls <- length(unique(unsampled$SPECIFICHAUL))
        
        cat("Year", i, "\tHauls", Nhauls, "\tUnsampled", Un.hauls, "\n"); flush.console();
        
        #Can only calculate a rate factor from the HAKE catch. Bycatch is the official total weight minus the hake catch (sp=206)
        #Use only sampled hauls and hake to calculate bycatch
        x <- x[(!is.na(x$HAUL_SAMPLED_BY)) & x$HAUL_SAMPLED_BY != 0 & x$SPECIES == species,]
        x$ByCatch <- x$OFFICIAL_TOTAL_CATCH - x$EXTRAPOLATED_WEIGHT 
        
        if (nrow(unsampled) > 0) { 
            # Use the sampled hauls to generate a monthly bycatch rate.
        
            B.byMonth <- aggregate(x$ByCatch, list(x$Month), sum, na.rm=T)
            T.byMonth <- aggregate(x$OFFICIAL_TOTAL_CATCH, list(x$Month), sum, na.rm=T)
            names(B.byMonth) <- c("Month", "ByCatch")
            names(T.byMonth) <- c("Month", "Catch")
            months <- unlist(T.byMonth[1])
            names(months) <- NULL
        
            bycatch.rate <- B.byMonth[,2] / T.byMonth[,2]
        
            Fraction.of.total <- rep(NA,12) #for each month
            Fraction.of.total[months] <- 1-bycatch.rate
        
            #Apply bycatch rate to unsampled hauls with total catch
            unsampled$ActualCatch <- unsampled$OFFICIAL_TOTAL_CATCH * Fraction.of.total[unsampled$Month] 
            unsampled$ByCatch <- unsampled$OFFICIAL_TOTAL_CATCH - unsampled$OFFICIAL_TOTAL_CATCH * Fraction.of.total[unsampled$Month] 
                
            # Put the dataset back together
            x <- rbind(x, unsampled)
        } # End if there are unsampled hauls.
        
        outdata <- data.frame(Sector="DomesticAtSea",Vessel=x$VESSEL, Month= x$Month, Day=x$Day, Year=i, lat=x$RETRV_LATITUDE,lon=x$RETRV_LONGITUDE,Catch.MT=x$ActualCatch/1000)
        #colnames(outdata) = c("Month", "Year", "Catch", "Sector")
        allOut <- rbind(allOut,outdata)
        
        if(!is.null(outfname)) {
            if (i == years[1]) {
                #    write.table(file=outfname, outdata, col.names=T, row.names=F, sep=",", append = T)
                write.csv(file=outfname, outdata, col.names=T, row.names=F, append = F)
            } else {
                write.csv(file=outfname, outdata, col.names=F, row.names=F, sep=",", append = T)
            } # End if-else
        }
    } # End for i
    
    allOut$lat <- convertToDecDeg(allOut$lat)
    allOut$lon <- -1*convertToDecDeg(allOut$lon)
    
    invisible(allOut)
}

convertToDecDeg <- function(x,rnd=2) {
#assumes that it is a four or five digit number, with the last two digits being minutes (0-59)
    deg <- floor(x/100)
    minutes <- 100*(x/100-deg)/60
    return(round(deg+minutes,rnd))
}




NorpacUnsampHauls <- function(ncatch,species=206) {
    ##########################
    #Return number of sampled and unsampled hauls by year
    #Code developed by Andi Stephens in 2010
    #Modified by Allan Hicks in 2013
    #
    #  Domestic Norpac catch consists of sampled hauls and unsampled hauls.
    #  For sampled hauls, EXTRAPOLATED_WEIGHT is the catch.  For unsampled
    #  hauls, we get a monthly bycatch rate and apply it to reduce the
    #  OFFICIAL_TOTAL_CATCH.
    #
    #  Assumptions:
    #
    #    All sampled hauls have an entry for SPECIES == 206, i.e., Hake.
    #    Each unsampled haul is represented by one unique record.
    ##############################
    
    colnames(ncatch)[1] <- "HAULJOIN"
    ncatch$SPECIFICHAUL <- paste(format(ncatch$HAULJOIN,digits=19),ncatch$HAUL,sep="_")
    ncatch$Month <- as.numeric(substring(ncatch$HAUL_DATE,6,7))
    ncatch$Year <- as.numeric(substring(ncatch$HAUL_DATE,1,4))

    years <- sort(unique(ncatch$Year))

    allOut <- NULL 
    for (i in years) {
        x <- ncatch[ncatch$Year == i,]  #only the specific year
        x <- x[(!is.na(x$OFFICIAL_TOTAL_CATCH)) & x$OFFICIAL_TOTAL_CATCH != 0,]  #only observations with an official catch
        
        Nrec <- nrow(x)
        Nhauls <- length(unique(x$SPECIFICHAUL))
        
        #Separate the unsampled hauls.  Only one record per haul if unsampled
        unsampled <- x[is.na(x$HAUL_SAMPLED_BY) | x$HAUL_SAMPLED_BY == 0,]
        Un.hauls <- length(unique(unsampled$SPECIFICHAUL))
        
        allOut <- rbind(allOut,data.frame(Year=i,Hauls=Nhauls, Unsampled=Un.hauls))
        
    } # End for i
    
    invisible(allOut)
}
