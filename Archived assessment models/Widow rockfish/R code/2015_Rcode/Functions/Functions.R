queryDB <- function(queryFilename,db,uid,sp="PWHT", start="2008", end="2009", querydir="sql/", asis=F) {

    require(RODBC)

    make.query = function(query, sp="PWHT", start="2008", end="2009", querydir="sql/") {
        #############################################################################
        #
        # Function make.query creates an SqlPlus query from a file.
        # Reads in query text and removes comment lines
        # Returns the query as one statement (list of characters) to pass to sqlQuery
        #
        # Andi Stephens, June 2010
        # reformatted and updated by Allan Hicks for Windows and RODBC, November 2011
        #
        #############################################################################
        # Get the query text
        qname <- paste(querydir,query,".query",sep="")
        if (!file.exists(qname)) { 
            stop("No such file: ", qname, "\n");
        }
        
        qq <- readLines(qname)
        
        # Make appropriate substitutions
        qq <- gsub("&sp", sp, qq)
        qq <- gsub("&beginyr", start, qq)
        qq <- gsub("&endyr", end, qq)
        
        #Remove comment lines
        qq <- qq[-grep("REM",qq)]  
        
        #Return query to be performed
        return(paste(qq,collapse=" "))
    }


    getData.fn <- function(dsn,uid,query,ais=asis) {
        ############################################################################    
        # Function to extract from a databse given a query
        ############################################################################
        
        require(RODBC)
        db <- odbcConnect(dsn=dsn,uid=uid)
        out <- sqlQuery(db, query,as.is=ais)
        odbcClose(db)
        
        return(out)
    }

    #make the query as a character string from the file in the sql directory
    query <- make.query(queryFilename,sp=sp,start=start,end=end,querydir=querydir)

    #get data from the database
    out <- getData.fn(dsn=db,uid=uid,query=query)

    return(out)
}


############################################################################### 
#
#   AUTHOR:  John R. Wallace (John.Wallace@noaa.gov) 
#   REVISED: Andi Stephens, 2010.
# 
#   Takes two tables with a shared primary key, and 
#   returns the rows of the second table for which the 
#   keys match those of the first.
#
#   NOTE:  The way this is used assumes that the second table is a
#          superset of the first (i.e., that each value is matched).
#
#   Changes:  
#
#        Changed name from original "match.f" to "find.matching.rows".
#
#        Removed sub-function 'paste.col' and made it standalone.
#
#        The matching function no longer modifies it's inputs, just
#        returns the values to be 'cbound' in the calling function.
#
##############################################################################


find.matching.rows <- function(file, table, findex = 1, tindex = 1, tcol = 2, round. = T) {
#############################################################################
# Function find.matching.rows
#
# Using the primary keys in columns named 'findex' and 'tindex', finds the
# matching values for 'file' in 'table' and returns 'table' column 'tcol'.
#
# Note that no test is made to determine if there are unmatched rows.
#############################################################################

  # Coerce a vector argument into a matrix
  if (is.null(dim(file))) {  dim(file) <- c(length(file), 1) }

  # If the primary keys are numeric, round them.
  if (round.) { 
    if (is.numeric(file[, findex])) { file[, findex] <- round(file[, findex]) }
    if (is.numeric(table[, tindex])) { table[, tindex] <- round(table[, tindex]) }
  } # End if round.

  # Convert the indices to character strings for comparison, and get the 
  # positions of the 'file' values in the 'table' values.
  matched.rows = match(paste.col(file[, findex]), paste.col(table[, tindex]))

  # Return the 'tcol' values in the rows of the 'table' that matched.
  return(table[matched.rows, tcol, drop = F]) 

} # End function find.matching.rows
 



paste.col <- function(x) {
#############################################################################
# Function paste.col
# 
# Converts a row to a string, "pasting" the columns together.
#############################################################################
  # If it's just a vector, return each value as a string.
  if (is.null(dim(x))) {
    return(paste(as.character(x))) 
  } # End if

  # Otherwise, it's a matrix.
  # Get the value of the first column in character form
  out <- paste(as.character(x[, 1])) 

  # Add on each of the rest of the columns
  for (i in 2:ncol(x)) { 
   out <- paste(out, as.character(x[, i])) 
  } # End for 

  return(out) 
} # End paste.col
 


try.detach = function(x) {
##############################################################################
#
# try.detach detaches iteratively and silently
#
##############################################################################
test = NULL
while (is.null(test)) {
  test = try(detach(deparse(substitute(x)), character.only=T), silent=T)
  } # End while
} # End function try.detach



replace.zeros = function(dest, source) {
##############################################################################
# Function replace.zeros replaces zero and NA values in one vector with
# values from a second.  Returns changes, and number of replacements;
# does NOT modify input.
# Function modified by ACH Nov 2011 to return the correct number of replacements when both vectors have NA's
##############################################################################
  # Replace NAs with zeros
  dest[is.na(dest)] = 0

  # Check for zeros and adjust
  test = (dest == 0)
  if (!is.na(sum(test))) {
    dest[test] = source[test]
  } # End if

  #There will be NA's if source has NA's where dest had zero or NA
  #test = sum(test) - sum(dest == 0)
  test = sum(test) - sum(is.na(dest))  #Sum of where zeros used to be, minus where NA's now are
  return(list(dest, sum(test)))
} # End function replace.zeros



extract.matrix = function(x) {
##############################################################################
# Function extract.matrix extracts a matrix from a table object.
##############################################################################

  return.val = NULL
  vector.names = NULL

  # Iterate over the rows of data,
  # binding them into a new matrix.
  tmp = attributes(x)$dim

  if (length(tmp) == 0) {
    # It's a vector, not a table
    return(x)
  } # End if

  if (length(tmp) > 1) {
    for ( i in 1:tmp[1] ) {
      return.val = rbind(return.val, x[i,])
    } # End for
    # Grab the row names.  The column names came along
    # automatically.
    rownames(return.val) = attributes(x)$dimnames[[1]]
  } else {
    # The table contains a vector of values (a one-D table)
    for ( i in 1:tmp[1] ) {
      return.val = c(return.val, x[[i]])
      vector.names = c(vector.names, attributes(x[i])$names)
    } # End for
    names(return.val) = vector.names
  } # End if-else

  return(return.val)
} # End extract.matrix


as.pctages = function(x) {
##############################################################################
#
# as.pctages converts the rows of a matrix to percentages cleanly.
#
##############################################################################
  # If it's a vector, convert to a matrix.
  if ( is.null(dim(x)) ) {
    dim(x) = c(1, length(x))
  } # End if

  return.val = NULL

  for ( i in 1:nrow(x) ) {
    if (sum(x[i,]) == 0) {
      tmp = rep(0, ncol(x)) 
    } else {
      tmp = x[i,]/sum(x[i,])
    } # End if-else
     return.val = rbind(return.val, tmp)
  } # End for i

  return(return.val)
} # End function as.pctages
