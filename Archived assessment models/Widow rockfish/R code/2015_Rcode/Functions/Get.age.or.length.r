##############################################################################
#
#  Get.age.or.length extracts age/length composition by weight from a table
#  assumed to be SEX x LENGTH x AGE and to determine the number of trips and
#  of fish associated with each length report from similarly shaped tables.
#
#  Behaviour controlled by 'output', which can be any of
#
#    age.and.length
#    age.only
#    length.only
#
##############################################################################

get.age.or.length = function(target, trips, fish, output, by.sex, pct) {

# Age-columns and length-rows

if ( output == "age.and.length") {

  n = attributes(trips)$dim[3]    # LENGTHs

  ftrips = rep(0, n)
  mtrips = rep(0, n)
  All.Trips = rep(0, n)

  ffish = rep(0, n)
  mfish = rep(0, n)
  All.Fish = rep(0,n)

  for ( i in 1:n ) {

    f_tdata = rowSums(trips[,1,i,])
    m_tdata = rowSums(trips[,2,i,])

    a_tdata = cbind(f_tdata,m_tdata)
    a_tdata[a_tdata > 1] = 1

    ftrips[i] = sum(a_tdata[,1])
    mtrips[i] = sum(a_tdata[,2])

    a_tdata = rowSums(a_tdata)
    a_tdata[a_tdata > 1] = 1

    All.Trips[i] = sum(a_tdata)

    ffish[i] = sum(fish[1,i,])
    mfish[i] = sum(fish[2,i,])
    All.Fish[i] = sum(fish[,i,])

  } # End for i

  if ( by.sex ) {

    fdata = extract.matrix(target[1,,])
    mdata = extract.matrix(target[2,,])

    my.colnames = colnames(fdata)
    f.rownames = rownames(fdata)
    m.rownames = rownames(mdata)

    NF = ncol(fdata)

    if (pct) {

      pct.data = round(as.pctages(cbind(fdata,mdata)), 6) 

      fdata = pct.data[,1:NF]
      mdata = pct.data[,(NF+1):(ncol(pct.data))]

    } # End if pct

    data = rbind(fdata, mdata)
    colnames(data) = my.colnames

    # Find and remove columns without value

    non.zero.rows = rowSums(data) > 0

    my.colnames = colnames(data)

    Gender = c(rep("Female", nrow(fdata)), rep("Male", nrow(mdata)))
    N.fish = c(ffish, mfish)
    N.trips = c(ftrips, mtrips)
    Len.bins = c(f.rownames, m.rownames)

    data = cbind(Gender, Len.bins, N.fish, N.trips, data)

    colnames(data) = c("Gender", "Len.bins", "N.fish", "N.trips", my.colnames)

    # Now remove rows without information

    data = data[non.zero.rows , ]

  } else {

    data = extract.matrix(target[1,,]) + extract.matrix(target[2,,])

    Len.bins = rownames(data)

    non.zero.rows = rowSums(data) > 0

    my.colnames = colnames(data)

    if (pct) { data = round(as.pctages(data), 6) }

    data = cbind(Len.bins, All.Fish, All.Trips, data)

    colnames(data) = c("Len.bins", "N.fish", "N.trips", my.colnames)

    data = data[non.zero.rows , ]

  } # End if-else by.sex

  return(data)

} # End age.and.length

# The length.only case

if ( output == "length.only" ) {

  n = attributes(trips)$dim[3]  # Lengths

  ftrips = NULL
  mtrips = NULL
  All.Trips = NULL

  ffish = rep(0, n)
  mfish = rep(0, n)
  All.Fish = rep(0,n)

  for ( i in 1:n ) {

    # Accumulate trip tables
    # These must be tables because
    # if we do the same thing as
    # above, we end up counting
    # multiple times.  Just believe
    # me on this one....

    if ( is.null(ftrips) ) {

       mtrips = trips[,2,i,] 
       ftrips = trips[,1,i,]
       All.Trips = trips[,,i,]

    } else {

      mtrips = mtrips + trips[,2,i,]  
      ftrips = ftrips + trips[,1,i,]
      All.Trips = All.Trips + trips[,,i,]

    } # End if-else

    ffish[i] = sum(fish[1,i,])
    mfish[i] = sum(fish[2,i,])
    All.Fish[i] = sum(fish[,i,])

  } # End for i

  mtrips = rowSums(mtrips)
  ftrips = rowSums(ftrips)
  All.Trips = rowSums(All.Trips)

  mtrips[mtrips > 1] = 1
  ftrips[ftrips > 1] = 1
  All.Trips[All.Trips > 1] = 1

  mtrips = sum(mtrips)
  ftrips = sum(ftrips)
  All.Trips = sum(All.Trips)
  
  ffish = sum(ffish)
  mfish = sum(mfish)
  All.Fish = sum(All.Fish)

  if ( by.sex ) {

    # Get age by length table of weights

    fdata = t(extract.matrix(target[1,,]))
    mdata = t(extract.matrix(target[2,,]))

    my.colnames = cbind(colnames(fdata), colnames(mdata))

    fdata = colSums(fdata)
    mdata = colSums(mdata)

    data = cbind(fdata,mdata)

    dim(data) = c(1,length(data))

    if (pct) { data = round(as.pctages(data), 6) }
    
    colnames(data) = my.colnames
    
    ntrips = All.Trips
    nfish = All.Fish

    data = cbind(nfish, ntrips, data)

    colnames(data) = c("N.fish","N.trips", my.colnames)

  } else {

    data = t(extract.matrix(target[1,,])) + t(extract.matrix(target[2,,]))

    my.colnames = colnames(data)

    data = colSums(data)

    if (pct) { data = round(as.pctages(data), 6) }

    ntrips = All.Trips
    nfish = All.Fish

    # Somehow necessary to get the colnames

    dim(data) = c(1,length(data))
    colnames(data) = my.colnames

    data = c(nfish, ntrips, data)
    dim(data) = c(1, length(data))
    colnames(data) = c("N.fish","N.trips", my.colnames)

  } # End if-else by.sex

  return(data)

} # End if length.only


# The age.only case

if ( output == "age.only" ) {

  n = attributes(trips)$dim[4]  # AGES
  ftrips = NULL
  mtrips = NULL
  All.Trips = NULL

  ffish = rep(0, n)
  mfish = rep(0, n)
  All.Fish = rep(0,n)

  for ( i in 1:n ) {

    # Accumulate trip tables

    if ( is.null(ftrips) ) {

       mtrips = trips[,2,,i] 
       ftrips = trips[,1,,i]
       All.Trips = trips[,,,i]

    } else {

      mtrips = mtrips + trips[,2,,i]  
      ftrips = ftrips + trips[,1,,i]
      All.Trips = All.Trips + trips[,,,i]

    } # End if-else

    ffish[i] = sum(fish[1,,i])
    mfish[i] = sum(fish[2,,i])
    All.Fish[i] = sum(fish[,,i])

  } # End for i


  mtrips = rowSums(mtrips)
  ftrips = rowSums(ftrips)
  All.Trips = rowSums(All.Trips)

  mtrips[mtrips > 1] = 1
  ftrips[ftrips > 1] = 1
  All.Trips[All.Trips > 1] = 1

  mtrips = sum(mtrips)
  ftrips = sum(ftrips)
  All.Trips = sum(All.Trips)
  
  ffish = sum(ffish)
  mfish = sum(mfish)
  All.Fish = sum(All.Fish)

  if ( by.sex ) {

    # Get age by length table of weights

    fdata = extract.matrix(target[1,,])
    mdata = extract.matrix(target[2,,])

    my.colnames = c(colnames(fdata),colnames(mdata))

    fdata = colSums(fdata)
    mdata = colSums(mdata)

    data = cbind(fdata,mdata)

    dim(data) = c(1,length(data))

    if (pct) { data = round(as.pctages(data), 6) }
    
    colnames(data) = my.colnames
    
    ntrips = All.Trips 
    nfish = All.Fish

    data = cbind(nfish, ntrips, data)

    colnames(data) = c("N.fish","N.trips", my.colnames)

  } else {

    data = extract.matrix(target[1,,]) + extract.matrix(target[2,,])

    my.colnames = colnames(data)

    data = colSums(data)

    if (pct) { data = round(as.pctages(data), 6) }

    ntrips = All.Trips
    nfish = All.Fish

    # Somehow necessary to get the colnames

    dim(data) = c(1,length(data))
    colnames(data) = my.colnames

    data = c(nfish, ntrips, data)
    dim(data) = c(1,length(data))
    colnames(data) = c("N.fish","N.trips", my.colnames)

  } # End if-else by.sex

  return(data)

} # End if age.only

} # End function get.age.or.length 

#
#  That's All, Folks!
#
##############################################################################
