SetUpWidowBDS.fn <- function(BDS,verbose=T,max.mmLength=1000,dataTypes=c("C"),
                                sampleMethods=c("R"),sampleTypes=c("","C","M"),
                                states=c("CA","OR","WA")) {
###############################################################################
#Reads in the BDS data and creates additional columns as well as filters the data
#Columns added are:
    #length.cm
#Start with:
#    load(PacFIN.REYE.bds.08.May.13.dmp)
#    BDS <- PacFIN.REYE.bds.08.May.13.dmp

###############################################################################


    BDS$length.cm <- BDS$FISH_LENGTH/10                                                 #lengths in cm

    BDS$inpfc <- rep(NA,nrow(BDS))                                  #assign single ipnfc codes to areas
    BDS$inpfc[BDS$INPFC_AREA %in% c("CA")] <- "CA"   #I put this in for CALCOM data
    BDS$inpfc[BDS$INPFC_AREA %in% c("CP")] <- "CP"
    BDS$inpfc[BDS$INPFC_AREA %in% c("MT")] <- "MT"
    BDS$inpfc[BDS$INPFC_AREA %in% c("EK","EU")] <- "EK"
    BDS$inpfc[BDS$INPFC_AREA %in% c("CL","COL","NC","SC")] <- "CL"
    BDS$inpfc[BDS$INPFC_AREA %in% c("VN","VUS")] <- "VN"
    ind1 <- BDS$inpfc %in% c("CA","CP","MT")  #South of 40.5
    ind2 <- BDS$inpfc %in% c("EK","CL","VN") #North of 40.5

    BDS$area  <- rep(NA,nrow(BDS))                                                      #set up area column for north and south, this will also classify out of stock areas
    BDS$area[ind1] <- "South"
    BDS$area[ind2] <- "North"
    if(verbose) {
        cat("Table of areas:")
        print(table(BDS$area))
        cat("\n")
        flush.console()
    }

    BDS$gear <- NA
    BDS$gear[BDS$GRID%in%c("SST","SHT","PWT","DST","DSG")] <- NA #"ShrimpTrawl"
    BDS$gear[BDS$GRID%in%c("RLT","GFT","GFS","GFL","FTS","FFT","BTT","BMT","56")] <- "BottomTrawl"
    BDS$gear[BDS$GRID%in%c("OTW","MDT","54")] <- "MidwaterTrawl"
    BDS$gear[BDS$GRID%in%c("PRT","DNT")] <- NA #"MiscTrawl"
    BDS$gear[BDS$GRID%in%c("BTR","CLP","CPT","FPT","OPT","PRW")] <- NA# "Pot"
    BDS$gear[BDS$GRID%in%c("HKL","JIG","LGL","OHL","POL","TRL","VHL")] <- "HnL"
    BDS$gear[BDS$GRID%in%c("DPN","DGN","GLN","ONT","SEN","STN")] <- "Net"
    BDS$gear[BDS$GRID%in%c("DVG","USP","MPT","UNK","XXX")] <- NA #"Other"   #MPT is CP midwater

    if(verbose) {
        cat("There are",sum(is.na(BDS$gear)),"rows out of",nrow(BDS),"where gear is not classified as",unique(BDS$gear),".\n\n")
        flush.console()
    }

    BDS$state <- rep(NA,nrow(BDS))
    BDS$state[BDS$SAMPLE_AGENCY == "CA"] <- "CA"
    BDS$state[BDS$SAMPLE_AGENCY == "OR"] <- "OR"
    BDS$state[BDS$SAMPLE_AGENCY == "W"] <- "WA"
    BDS$state[BDS$SAMPLE_AGENCY == "PW"] <- "PW"


    #drop columns that have all NA's
    ind <- apply(BDS,2,function(x){all(is.na(x))})  #goes column by column and returns TRUE for columns with all NA's
    if(verbose) {
        cat("These columns were dropped:\n")
        print(names(BDS)[ind])
        cat("\n")
        flush.console()
    }
    BDS <- BDS[,!ind]

    if(verbose) {
        cat("\nThe dimensions of the unfiltered dataframe are:\n")
        print(dim(BDS))
        cat("\n")
        flush.console()
    }

    #filter out the rows that will not be used
    keep <- rep(T,nrow(BDS))                                             #TRUE means to omit that row
    totKeep <- sum(keep)
    keep <- keep & !is.na(BDS$area)
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue of NA in area (Canadian catch)\n")}; totKeep <- sum(keep)
    keep <- keep & (BDS$state %in% states)
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue not a specified state:",states,"\n")}; totKeep <- sum(keep)
    keep <- keep & (!is.na(BDS$gear))
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue not a specified gear:",unique(BDS$gear),"\n")}; totKeep <- sum(keep)
    keep <- keep & !is.na(BDS$FISH_LENGTH)                               #omit missing observations of length
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue of NA in Fish Length\n")}; totKeep <- sum(keep)
    keep <- keep & BDS$FISH_LENGTH > 0                                   #omit lengths 0 or less
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue of Fish Length <= 0\n")}; totKeep <- sum(keep)
    keep <- keep & BDS$FISH_LENGTH <= max.mmLength                       #omit fish greater than maxLength
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue of Fish Length >",max.mmLength,"mm\n")}; totKeep <- sum(keep)
    keep <- keep & BDS$DATA_TYPE %in% dataTypes                          #types of data such as commercial or survey
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue dataTypes were not of type",dataTypes,"\n")}; totKeep <- sum(keep)

    midFilterBDS <- BDS[keep,]                                                  #save this to use for calcualting the median landing weights

    keep <- keep & BDS$SAMPLE_METHOD %in% sampleMethods                  #sampling method such as random, stratified, special,...
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue sampleMethods were not of type",sampleMethods,"\n")}; totKeep <- sum(keep)
    keep <- keep & BDS$SAMPLE_TYPE %in% sampleTypes                      #sampling type such as market, research, special
     if(verbose){cat(totKeep-sum(keep),"rows omitted becasue sampleTypes were not of type",sampleTypes,"\n")}; totKeep <- sum(keep)
     if(verbose) {
        cat(sum(keep),"rows retained out of",nrow(BDS),"based on filtering.\n")
     }
    BDS <- BDS[keep,]

    #BDS$length.cm <- round(BDS$length.cm,0)                         #round all lengths to nearest cm
    #if(verbose) cat("Lengths were rounded to the nearest cm in length.cm\n")
    BDS$length.cm <- floor(BDS$length.cm)                         #floor all lengths to nearest cm
    if(verbose) cat("Lengths were floored to the lower cm in length.cm\n")

    #set up total weights for expansion
    #use EXP_WT if provided (for OR only)
    #use TOTAL_WGT otherwise
    BDS$totalWt <- rep(NA,nrow(BDS))
    ind <- !is.na(BDS$TOTAL_WGT) & BDS$TOTAL_WGT>0                                      #use total weight  where available
    BDS$totalWt[ind] <- BDS$TOTAL_WGT[ind]
    ind <- !is.na(BDS$EXP_WT) & BDS$EXP_WT>0                                           #OR samples with an expanded total weight
    BDS$totalWt[ind] <- BDS$EXP_WT[ind]                                          # replace any total weights with these expanded (IAN used exp_wts for Ebglish Sole, but did not use TOTAL_WGT from OR)
    if(verbose) cat("A column call totalWt created using TOTAL_WGT and EXP_WT\n")

    if(verbose) {
        cat("\nThe final dimensions of the filtered dataframe are:\n")
        print(dim(BDS))
        cat("\n")
        flush.console()
    }

    BDS$SAMPLE_NO <- as.character(BDS$SAMPLE_NO)

    return(BDS)
}





# match.f <- function(file, table, findex = 1, tindex = findex, tcol = 2, round. = T)
# {
# #
# #   DATE WRITTEN:  ???      LAST REVISED:   7 March 2000
# #   AUTHOR:  John R. Wallace (John.Wallace@noaa.gov)
# #
#        paste.col <- function(x)
#        {
#                if(is.null(dim(x)))
#                        return(paste(as.character(x)))
#                out <- paste(as.character(x[, 1]))
#                for(i in 2:ncol(x)) {
#                        out <- paste(out, as.character(x[, i]))
#                }
#                out
#        }
#        if(is.null(dim(file))) {
#                dim(file) <- c(length(file), 1)
#        }
#        if(round.) {
#                for(i in findex) {
#                        if(is.numeric(file[, i]))
#                                file[, i] <- round(file[, i])
#                }
#                for(i in tindex) {
#                        if(is.numeric(table[, i]))
#                                table[, i] <- round(table[, i])
#                }
#        }
#        cbind(file, table[match(paste.col(file[, findex]), paste.col(table[, tindex])), tcol, drop = F])
# }
