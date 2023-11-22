# TODO: List of parameters that are given by default
# These might not go well for use with evens that have another timescale
# For example: behavior
# 
# ISITh = 100 msec
# cuttof = 0.1 sec
# These lines
# isi <- diff(st) * 1000
# cutoff <- cutoff * 1000
# eps <- 10^-10
# The params to merge the bursts given to the smejea function

# TODO: nomenclature
# I believe we should keep the "spike" as the main event
# if using for something else, it's easier to adapt than to
# change all definitions (e.g., inter-spike interval will become inter-event interval, lame)

# helper functions --------------------------------------------------------
#Finds peaks in logISI histogram
# The location of the peaks of this histogram are found using a custom peak finding algorithm
# V. Pasquale, S. Martinoia, and M. Chiappalone. A self-
# adapting approach for the detection of bursts and network
# bursts in neuronal cultures. J. Comput. Neurosci., 29(1-2):213–
# 229, 2010.
# TODO: pracma::findpeaks() performs slightly differently 
# Finding peaks is a bit noisy in general
get_peaks <- function(h,
                      Pd = 2,
                      Th = 0,
                      Np = NULL) {
  m <- 0
  # This is expecting a histogram object
  # with h$density that has been modified as h$counts/sum(h$counts)
  L <- length(h$density)
  j <- 0
  Np <- ifelse(is.null(Np), L, Np)
  pks <- NULL
  locs <- NULL
  void.th <- 0.7
  while ((j < L) && (m < Np)) {
    j <- j + 1
    endL <- max(1, j - Pd)
    if (m > 0 && j < min(c(locs[m] + Pd, L - 1))) {
      j <- min(c(locs[m] + Pd, L - 1))
      endL <- j - Pd
    }
    endR <- min(L, j + Pd)
    temp <- h$density[endL:endR]
    aa <- which(j == endL:endR)
    temp[aa] <- -Inf
    if (Pd > 1) {
      idx1 <- max(1, aa - 2)
      idx2 <- min(aa + 2, length(temp))
      idx3 <- max(1, aa - 1)
      idx4 <- min(aa + 1, length(temp))
      if (sum((h$density[j] > (temp[c(1:idx1, idx2:length(temp))] + Th)) ==
              FALSE) == 0 &&
          sum((h$density[j] > (temp[idx3:idx4])) == FALSE) == 0 &&
          j != 1 && j != L) {
        m <- m + 1
        pks[m] <- h$density[j]
        locs[m] <- j
      }
    } else if (sum((h$density[j] > (temp + Th)) == FALSE) == 0) {
      m <- m + 1
      pks[m] <- h$density[j]
      locs[m] <- j
    }
    
  }
  ret <- data.frame(pks = pks, locs = locs)
  return(ret)
}

#Function to find cutoff threshold.
find_thresh <- function(h, ISITh = 100) {
  # From the method definition:
  # The smallest ISI min i for which void(i)>0.7 
  # is set as the threshold for the maximum ISI 
  # in a burst
  void_th <- 0.7
  gp <- get_peaks(h)
  num_peaks <- nrow(gp)
  # The two lines below will fail if we don't check this
  stopifnot(num_peaks > 1)
  pkx <- h$breaks[gp$locs]
  # If no peak is found in the histogram with ISI ≤ ISIth, 
  # the spike train is classified as containing no bursts.
  # TODO: <= was not implemented in original code
  intra_indx <- which(pkx < ISITh)
  # which might return integer(0)
  # use this if/else to deal with that
  if (length(intra_indx) >= 1) {
    intra_pks <- gp$pks[intra_indx]
    # get the closes peak to the break
    max_intra <- max(intra_pks)
    # if points are sorted which.max() should be same as length()
    # which.max handles ties differently than length()
    max_idx <- which.max(intra_pks)
  } else {
    # if no pkx < ISIth return -1000
    return(-1000)
  }
  # get the values for the point of break
  x1 <- pkx[max_idx]
  # The largest peak of the histogram corresponding to an ISI less than or
  # equal to MCV is set as the intraburst peak, CIBP
  # CIBP = y1
  y1 <- max_intra
  locs1 <- gp$locs[max_idx]
  num_peaks_after_burst <- num_peaks - max_idx
  if (num_peaks_after_burst == 0) {
    return(NULL)
  } else {
    # subset detected peaks after the breaking point
    gp2 <- gp[(max_idx + 1):num_peaks, ]
    # We will iterate over these points because
    # We are looking for
    # the minimum value of the histogram between the 
    # intraburst peak (CIBP) and each of the fol-
    # lowing peaks, Cp i (i = 1, ..., N)
    ymin <- sapply(gp2$locs, function(xx)
      min(h$density[locs1:xx]))
    xmin <-
      sapply(gp2$locs, function(xx)
        which.min(h$density[locs1:xx])) + locs1 - 1
    # The method's definition
    # voidParameter = 1 - (Cmin)/(sqrt(CIBP)*Cpi)
    voidParameter <- 1 - (ymin / sqrt(y1 * gp2$pks))
  }
  # Find the first point that meets void_th criteria
  indxvoid <- suppressWarnings(min(which(voidParameter >= void_th)))
  if (is.infinite(indxvoid)) {
    flags <- c(1, 0)
    return(NULL)
  } else {
    # find the break point 
    # This point will separate intra and inter burst events
    break_point <- xmin[indxvoid]
    # find it in terms of h breaks
    ISImax <- h$breaks[break_point]
    return(ISImax)
  }
}

# Calculates value for burst detection by
# finding a threshold for cutting the histogram
# this will define itraburst and interburst intervals
logisi_break_calc <- function(st, cutoff) {
  # TODO: check where this 1000 is coming from
  # guess is that is a sec to ms conversion
  isi <- diff(st) * 1000
  cutoff <- cutoff * 1000
  max_isi <- ceiling(log10(max(isi)))
  # same here, greater than 1 is contaminated 
  # by the * 1000 
  # not including all events, only those with enough isi
  # if this is not subseted, histogram call below returns
  # Error in hist.default(isi, breaks = br, plot = FALSE) : 
  #  some 'x' not counted; maybe 'breaks' do not span range of 'x'
  isi <- isi[isi >= 1]
  # create logspace breaks for counting
  # not sure why we want it to go to 10*max_isi
  br <- pracma::logspace(0, max_isi, 10 * max_isi)
  h <- hist(isi, breaks = br, plot = FALSE)
  h$density <- h$counts / sum(h$counts)
  # smooth the density with 0.05 span, subset predicted y
  # this is quite overfitting
  h$density <- lowess(h$density, f = 0.05)$y
  thr <- find_thresh(h, cutoff)
  if (!is.null(thr)) {
    # scale back 
    thr <- thr / 1000
  }
  thr
}

###Function to add burst related spikes to edges of bursts
add.brs <- function(bursts, brs, spike_train) {
  burst.adj <-
    data.frame(beg = rep(0, dim(bursts)[1]), end = rep(0, dim(bursts)[1]))
  for (i in 1:dim(bursts)[1]) {
    for (j in 1:dim(brs)[1]) {
      if (data.table::between(bursts[i, 1], brs[j, 1], brs[j, 2]) |
          data.table::between(bursts[i, 2], brs[j, 1], brs[j, 2]))
      {
        burst.adj$beg[i] <- min(bursts[i, 1], brs[j, 1])
        burst.adj$end[i] <- max(bursts[i, 2], brs[j, 2])
        break
      } else {
        burst.adj$beg[i] <- bursts[i, 1]
        burst.adj$end[i] <- bursts[i, 2]
      }
      if (brs[j, 2] > bursts[i, 2]) {
        break
      }
    }
  }
  
  diff.begs <- diff(burst.adj[, "beg"])
  rep.bursts.begs <- which(diff.begs == 0)
  if (any(rep.bursts.begs)) {
    burst.adj <- burst.adj[-rep.bursts.begs, ]
  }
  diff.ends <- diff(burst.adj[, "end"])
  rep.bursts.end <- which(diff.ends == 0) + 1
  if (any(rep.bursts.end)) {
    burst.adj <- burst.adj[-rep.bursts.end, ]
  }
  start.times <- spike_train[burst.adj$beg]
  end.times <- spike_train[burst.adj$end]
  durn <- end.times - start.times
  len <- burst.adj$end - burst.adj$beg + 1
  mean.isis <- durn / (len - 1)
  N.burst <- dim(burst.adj)[1]
  IBI <- c(NA, start.times[-1] - end.times[-N.burst])
  result <-
    cbind(
      beg = burst.adj$beg,
      end = burst.adj$end,
      IBI = IBI,
      len = len,
      durn = durn,
      mean.isis = mean.isis,
      SI = rep(1, N.burst)
    )
  result
  
}

# These functions come from sjemea package
# It's somewhere in this generic file
# https://github.com/sje30/sjemea/blob/master/R/surprise2.R
num_bursts <- function(b) {
  ## Return the number of bursts found for a spike train.
  if(is.na(b[1])) {
    return(0)
  }  else {
    return(nrow(b))
  }
}

calc_ibi <- function(spikes, bursts) {
  ## Compute the interburst intervals (IBI) for one spike train.
  ## Only valid if more than one burst.
  nburst <- num_bursts(bursts)
  if ( nburst == 0) {
    #no bursts
    res <- NA                            
  } else {
    if (nburst == 1) {
      #cannot compute  IBI w/only 1 burst.
      res <- NA                          
    } else {
      ## find end spike in each burst.
      end <- bursts[,"beg"] + bursts[,"len"] - 1
      ## for NBURST bursts, there will be NBURST-1 IBIs.
      start_spikes <- bursts[2:nburst,"beg"]
      end_spikes <- end[1:(nburst-1)]
      ## NEX uses a strange definition of IBI -- it counts the time between
      ## the first spike of burst N and the first spike of burst N+1 as the
      ## IBI.  If we want to use that definition, use the following line:
      ##end.spikes   = b[1:(nburst-1),"beg"]
      res <- spikes[start_spikes] - spikes[end_spikes]
    }
  }
  return(res)
}

##Function for finding bursts, taken from sjemea
logisi_find_burst <- function(spikes, params, verbose = FALSE) {
  ## For one spike train, find the burst using log isi method.
  ## e.g.
  ## find.bursts(s$spikes[[5]])
  ## init. params currently in logisi_par
  ##
  ##beg.isi =    params$beg.isi
  ##end.isi =    params$end.isi
  #min_ibi <- params$min_ibi
  #min_durn <- params$min_durn
  #min_spikes <- params$min_spikes
  #isi_low <- params$isi_low
  
  #value to return if no bursts found.
  no_bursts <- NA
  # This will make the variables available in f(x) environment
  # We have to be somewhat cautious about the variables we pass 
  # to params list
  list2env(params, envir = environment())
  nspikes <- length(spikes)
  if (verbose) {
    usethis::ui_info("Running `logisi_find_burst` with params")
    usethis::ui_line("{names(params)}:{params}")
    usethis::ui_line("Number of events is: {usethis::ui_value(nspikes)}")
  }
  ## Create a temp array for the storage of the bursts.  Assume that
  ## it will not be longer than Nspikes/2 since we need at least two
  ## spikes to be in a burst.
  max_bursts <- floor(nspikes / 2)
  bursts <- matrix(NA, nrow = max_bursts, ncol = 3)
  colnames(bursts) = c("beg", "end", "IBI")
  # initialize burst number
  burst <- 0                            
  
  ## Phase 1 -- burst detection. ------------
  ## Each interspike interval of the data
  ## is compared with the threshold THRE. 
  # If the interval > THRE, it can not be part of a burst; 
  # else the interval may be part of a burst.
  
  # init variables
  ## last_end is the time of the last spike in the previous burst.
  ## This is used to calculate the IBI.
  last_end <- NA
  # eps is basically zero ?
  eps <- 10 ^ (-10)
  in_burst <- FALSE
  #for first burst, there is no IBI, start at n=2
  n <- 2
  # spike assignment into bursts
  # this process is iterative/cumulative, so we need a loop
  while (n < nspikes) {
    next_isi <- spikes[n] - spikes[n - 1]
    if (in_burst) {
      if (next_isi - isi_low > eps) {
        ## end of burst
        end <- n - 1
        in_burst <- FALSE
        ibi <-  spikes[beg] - last_end
        last_end <- spikes[end]
        res <- c(beg, end, ibi)
        # increase burst counter
        burst <- burst + 1
        # TODO: I think this has to be eliminated
        # browser() will give you nothing useful (?)
        if (burst > max_bursts) {
          print("too many bursts!!!")
          browser()
        }
        # append the results to the bursts matrix
        bursts[burst, ] <- res
      }
    } else {
      ## not yet in burst.
      if (next_isi - isi_low <= eps) {
        ## Found the start of a new burst.
        beg = n - 1
        in_burst = TRUE
      }
    }
    n <- n + 1
  }
  
  ## At the end of the burst, check if we were in a burst when the
  ## train finished.
  if (in_burst) {
    end <- nspikes
    ibi <-  spikes[beg] - last_end
    res <- c(beg, end, ibi)
    burst <- burst + 1
    if (burst > max.bursts) {
      print("too many bursts!!!")
      browser()
    }
    bursts[burst, ] <- res
  }
  
  ## Check if any bursts were found.
  if (burst > 0) {
    ## truncate to right length, 
    # as bursts will typically be very long.
    bursts = bursts[1:burst, , drop = FALSE]
  } else {
    ## no bursts were found, so return an empty structure.
    return(no_bursts)
  }
  
  if (verbose) {
    usethis::ui_done("End of phase1")
    usethis::ui_line("{usethis::ui_value(burst)} bursts detected")
    usethis::ui_info("Start of Phase2: Merging bursts")
  }
  
  
  ## Phase 2 -- merging of bursts.  ------
  # Here we see if any pair of bursts have an IBI < MIN.IBI;
  # if so, we then merge the bursts.
  ## We specifically need to check when say three bursts 
  # are merged into one.
  # TODO: there might be a tidyverse way of doing this
  # but don't fix what's not broken. 
  # Also, bursts is a matrix not a data.frame
  
  ibis <- bursts[, "IBI"]
  merge_bursts <- which(ibis < min_ibi)
  
  if (any(merge_bursts)) {
    ## Merge bursts efficiently.  
    # Work backwards through the list, and
    ## then delete the merged lines afterwards.  
    # This works when we
    ## have say 3+ consecutive bursts that merge into one.
    
    for (burst in rev(merge_bursts)) {
      # assign the previous end value to the burst value
      bursts[burst - 1, "end"] <- bursts[burst, "end"]
      # not needed, but helpful. 
      # could be used to subset by doing complete.cases()
      bursts[burst, "end"] <- NA         
    }
    # remove the unwanted info.
    bursts <- bursts[-merge_bursts, , drop = FALSE] 
    if (verbose) {
      usethis::ui_line("{usethis::ui_value(length(merge_bursts))} bursts were merged")
    }
  }
  
  if (verbose) {
    usethis::ui_done("End of Phase 2")
    usethis::ui_info("Start of Phase 3: remove small bursts")
  }
  
  
  ## Phase 3 -- remove small bursts: less than min duration 
  ## (MIN.DURN), or having too few spikes (less than MIN.SPIKES).
  ## CAUTION: In this phase,
  ## we have the possibility of deleting all spikes !!!!
  
  ## LEN = number of spikes in a burst.
  ## DURN = duration of burst.
  len <- bursts[, "end"] - bursts[, "beg"] + 1
  durn <- spikes[bursts[, "end"]] - spikes[bursts[, "beg"]]
  # add new data to bursts matrix
  bursts <- cbind(bursts, len, durn)
  ## we have the possibility of deleting all spikes !!!!
  rejects <- which ((durn < min_durn) | (len < min_spikes))
  if (any(rejects)) {
    if (verbose) {
      usethis::ui_warn("{usethis::ui_value(length(rejects))} out of {nrow(bursts)} bursts removed because criteria not met")
    }
    bursts <- bursts[-rejects, , drop = FALSE]
  }
  
  if (nrow(bursts) == 0) {
    ## All the bursts were removed during phase 3.
    bursts <- no_bursts
  } else {
    ## Compute mean ISIS
    len <- bursts[, "end"] - bursts[, "beg"] + 1
    durn <- spikes[bursts[, "end"]] - spikes[bursts[, "beg"]]
    mean.isis <- durn / (len - 1)
    
    ## Recompute IBI (only needed if phase 3 deleted some cells).
    if (nrow(bursts) > 1) {
      ibi2 <- c(NA, calc_ibi(spikes, bursts))
    } else {
      ibi2 <- NA
    }
    bursts[, "IBI"] <- ibi2
    
    SI <- rep(1, length(mean.isis))
    bursts <- cbind(bursts, mean.isis, SI)
  }
  
  ## End -- return burst structure. ----
  if (verbose) {
    usethis::ui_done("End of Phase 3")
  }
  return(bursts)
}


# logISI pasq method ------
logisi_pasq_method <- function(spike_train, cutoff = 0.1) {
  cutoff <- ifelse(is.null(cutoff), 0.1, cutoff)
  if (length(spike_train) > 3) {
    #Calculates threshold as isi_low
    isi_low <- logisi_break_calc(spike_train, cutoff)
    if (is.null(isi_low) || isi_low >= 1) {
      #If no value for isi_low found, 
      # or isi_low above 1 second, 
      # find bursts using threshold equal to cutoff 
      # (default 100ms)
      logisi_par <- list(
        min_ibi = 0,
        min_durn = 0,
        min_spikes = 3,
        isi_low = cutoff
      ) 
      result <- logisi_find_burst(spike_train, logisi_par)
    } else if (isi_low < 0) {
      result <- NA
    } else if (isi_low > cutoff & isi_low < 1) {
      logisi_par <- list(
        min.ibi = isi_low,
        min.durn = 0,
        min.spikes = 3,
        isi_low = cutoff
      ) #If isi_low >cutoff, find bursts using threshold equal to cutoff (default 100ms)
      bursts <- logisi_find_burst(spike_train, logisi_par)
      if (!is.na(bursts)[1]) {
        logisi_par2 <- list(
          min.ibi = 0,
          min.durn = 0,
          min.spikes = 3,
          isi_low = isi_low
        ) #If bursts have been found, add burst related spikes using threshold of isi_low
        brs <- logisi_find_burst(spike_train, logisi_par2)
        result <- add.brs(bursts, brs, spike_train)
      } else {
        result <- bursts
      }
    } else {
      logisi_par <- list(
        min.ibi = 0,
        min.durn = 0,
        min.spikes = 3,
        isi_low = isi_low
      ) #If isi_low<cutoff, find bursts using a threshold equal to isi_low
      result <- logisi_find_burst(spike_train, logisi_par)
    }
    
  } else {
    result <- NA
  }
  result
}

