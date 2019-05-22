library(microbenchmark)
library(parallel)

# there was a bug in the original R function from the seewave package, so I copied it here and fixed it
# (there was an overflow, probably because it wasn't made to use large time-bins for the fft)
myspec <-
  function (wave,
            f,
            wl = 512,
            wn = "hanning",
            fftw = FALSE,
            norm = TRUE,
            scaled = FALSE,
            PSD = FALSE,
            PMF = FALSE,
            correction = "none",
            dB = NULL,
            dBref = NULL,
            at = NULL,
            from = NULL,
            to = NULL,
            identify = FALSE,
            col = "black",
            cex = 1,
            plot = 1,
            flab = "Frequency (kHz)",
            alab = "Amplitude",
            flim = NULL,
            alim = NULL,
            type = "l",
            ...)
  {
    if (!isTRUE(norm) & PMF)
      stop("'PMF' can be computed only if 'norm' is TRUE")
    if (!isTRUE(norm) & !is.null(dB))
      stop("dB are computed on normalised spectra only, 'norm' should be turned to TRUE")
    if (norm & scaled)
      stop("'norm' and 'scaled' cannot be both set to TRUE")
    if (PMF & scaled)
      stop("'norm' and 'PMF' cannot be both set to TRUE")
    if (!is.null(dB) & PMF)
      stop("PMF cannot be in dB")
    if (is.null(dB) & !is.null(dBref))
      stop("'dB' cannot be NULL  when 'dBref' is not NULL")
    if (is.logical(dB))
      stop("'dB' is no more a logical. Please see the documentation: help(spec).")
    if (!is.null(dB) && all(dB != c("max0", "A", "B", "C", "D")))
      stop("'dB' has to be one of the following character strings: 'max0', 'A', 'B', 'C' or 'D'")
    input <- inputw(wave = wave, f = f)
    wave <- input$w
    f <- input$f
    rm(input)
    if (!is.null(from) | !is.null(to)) {
      if (is.null(from) && !is.null(to)) {
        a <- 1
        b <- round(to * f)
      }
      if (!is.null(from) && is.null(to)) {
        a <- round(from * f)
        b <- length(wave)
      }
      if (!is.null(from) && !is.null(to)) {
        if (from > to)
          stop("'from' cannot be superior to 'to'")
        if (from == 0) {
          a <- 1
        }
        else
          a <- round(from * f)
        b <- round(to * f)
      }
      wl <- (b - a) + 1
      wave <- as.matrix(wave[a:b, ])
    }
    if (!is.null(at)) {
      c <- round(at * f)
      wl2 <- wl %/% 2
      wave <- as.matrix(wave[(c - wl2):(c + wl2), ])
    }
    n <- nrow(wave)
    W <- ftwindow(n, wn = wn, correction = correction)
    wave <- wave * W
    if (fftw == FALSE) {
      y <- Mod(fft(wave[, 1]))
    }
    else {
      p <- fftw::planFFT(n)
      y <- Mod(fftw::FFT(wave[, 1], plan = p))
    }
    if (scaled) {
      y <- y / length(y)
    }
    y <- 2 * y[1:(n %/% 2)]
    if (PSD) {
      y <- y ^ 2
    }
    if (norm) {
      y <- y / max(y)
    }
    if (PMF) {
      y <- y / sum(y)
    }
    if (!is.null(dB)) {
      y <- ifelse(y == 0, yes = 1e-06, no = y)
      if (is.null(dBref)) {
        y <- 20 * log10(y)
      }
      else {
        y <- 20 * log10(y / dBref)
      }
      if (dB != "max0") {
        if (dB == "A")
          y <- dBweight(x * 1000, dBref = y)$A
        if (dB == "B")
          y <- dBweight(x * 1000, dBref = y)$B
        if (dB == "C")
          y <- dBweight(x * 1000, dBref = y)$C
        if (dB == "D")
          y <- dBweight(x * 1000, dBref = y)$D
      }
    }
    if (is.null(alim)) {
      if (is.null(dB)) {
        alim <- c(0, 1.1)
      }
      else {
        alim <- c(min(y, na.rm = TRUE), max(y, na.rm = TRUE) +
                    20)
      }
      if (PMF | !isTRUE(norm))
        alim <- c(0, max(y, na.rm = TRUE))
    }
    if (is.null(at)) {
      multiplicator <- f / n / 1000 # added this line to prevent overflow
      x <- ((0:(n - 1)) * multiplicator)[1:(n %/% 2)]
    }
    else {
      x <- ((0:(wl - 1)) * f / wl / 1000)[1:(wl %/% 2)]
    }
    if (plot == 1) {
      if (!is.null(dB)) {
        plot(
          x = x,
          y = y,
          xaxs = "i",
          xlab = flab,
          xlim = flim,
          yaxs = "i",
          yaxt = "s",
          ylab = alab,
          ylim = alim,
          col = col,
          cex = cex,
          type = type,
          las = 1,
          ...
        )
      }
      else {
        if (isTRUE(norm)) {
          yaxt <- "n"
          ylab <- alab
          if (isTRUE(PMF)) {
            yaxt = "s"
          }
        }
        else {
          yaxt <- "s"
          ylab <- " "
        }
        plot(
          x = x,
          y = y,
          xaxs = "i",
          xlab = flab,
          xlim = flim,
          yaxs = "i",
          yaxt = yaxt,
          ylab = ylab,
          ylim = alim,
          col = col,
          cex = cex,
          type = type,
          las = 1,
          ...
        )
      }
      if (identify) {
        cat("Choose points on the spectrum\n")
        if (.Platform$OS.type == "windows")
          flush.console()
        id <- identify(
          x = x,
          y = y,
          labels = round(x, 2),
          tolerance = 0.15,
          col = "red"
        )
        id.freq <- x[id]
        id.amp <- y[id]
        coord <- list(freq = id.freq, amp = id.amp)
        return(coord)
      }
    }
    if (plot == 2) {
      if (!is.null(dB)) {
        plot(
          x = y,
          y = x,
          xaxs = "i",
          xlab = alab,
          xlim = alim,
          yaxs = "i",
          yaxt = "s",
          ylab = flab,
          ylim = flim,
          col = col,
          cex = cex,
          type = type,
          las = 1,
          ...
        )
      }
      else {
        if (isTRUE(norm)) {
          xaxt <- "n"
          xlab <- alab
          if (isTRUE(PMF)) {
            xaxt = "s"
          }
        }
        else {
          xaxt <- "s"
          xlab <- " "
        }
        plot(
          x = y,
          y = x,
          xaxs = "i",
          xaxt = xaxt,
          xlab = xlab,
          xlim = alim,
          yaxs = "i",
          ylab = flab,
          ylim = flim,
          col = col,
          cex = cex,
          type = type,
          las = 1,
          ...
        )
      }
      if (identify) {
        cat("choose points on the spectrum\n")
        if (.Platform$OS.type == "windows")
          flush.console()
        id <- identify(
          x = y,
          y = x,
          labels = round(x, 2),
          tolerance = 0.15,
          col = "red"
        )
        id.freq <- x[id]
        id.amp <- y[id]
        coord <- list(freq = id.freq, amp = id.amp)
        return(coord)
      }
    }
    if (plot == 1 | plot == 2) {
      spec <- cbind(x, y)
      invisible(spec)
    }
    else if (plot == FALSE) {
      spec <- cbind(x, y)
      return(spec)
    }
  }

substrRight <- function(x, n) {
  # helper function to get the substring from the rigth
  # Params:
  #  x: string
  #  n: number of chars to take
  # Returns:
  #  the substring
  substr(x, nchar(x) - n + 1, nchar(x))
}


get.timestamp <- function(file) {
  # extract the timestamp from the filename
  # in case of the be a bee data the filenames are in a form like
  # "miks boden-03-may-20-52.wav"
  # as I was using a script to convert mp3 back to wav, some files end with
  # mp3.wav
  # Params:
  #  file: the file-name
  # Returns: the timestamp
  last <- substrRight(file, 16)
  # in case of convertion
  if (substrRight(last, 7) == "mp3.wav") {
    last <- substrRight(file, 20)
  }
  day <- substr(last, 1, 2)
  month.str <- substr(last, 4, 6)
  if (month.str == "jun") {
    month = "06"
  } else if (month.str == "may") {
    month = "05"
  } else if (month.str == "apr") {
    month = "04"
  }
  year <- "2013"
  hour <- substr(last, 8, 9)
  minute <- substr(last, 11, 12)
  seconds <- "00"
  date <- paste(year, month, day, sep = "-")
  time <- paste(hour, minute, sep = ":")
  return(as.POSIXlt(paste(date, time)))
}

quoted <- function(s) {
  return(paste("\"", s, "\"", sep = ""))
}

analyse <-
  function(file,
           directory,
           seconds,
           channel.number,
           csv,
           maxfreq = 1000,
           binfreq = 20) {
    # Run a fft on the file and get some addtional information like dominant frequency.
    # Saves the extracted information to specified csv file
    #
    # Params:
    #  file: the file name
    #  directory: the path to the file
    #  seconds: bin width for the fft
    #  channel.number: 1 for left channel, 2 for rigth channel
    #  csv: the name of an existing csv-file to save the results to
    #  maxfreq: maximal frequency to save to the csv
    #  binfreq: width of the frequency bands
    
    # read header of wav only, to obtain duration
    
    print(file)
    wav <- readWave(paste(directory, file, sep = ""), header = TRUE)
    samp.rate <- wav$sample.rate
    duration <- round(wav$samples / samp.rate, 0)
    # Devide in Chunks
    # run only if long enough
    if (duration < seconds) {
      stop("File to short")
    }
    if (channel.number < 0 || channel.number > 2) {
      stop("Channel number must be 1 or 2")
    }
    if (channel.number > wav$channels) {
      stop("Channel does not exist")
    }
    froms <- seq(0, duration - seconds, seconds)
    current.timestamp <- get.timestamp(file)
    rest <- minute(current.timestamp) %% 10
    if (rest != 0) {
      froms <- froms + 60 * (10 - rest)
      current.timestamp <- current.timestamp + 60 * (10 - rest)
    }
    for (i in 1:length(froms)) {
      to <- froms[i] + seconds
      if (to > duration) {
        to <- duration
      }
      wav <-
        readWave(
          paste(directory, file, sep = ""),
          from = froms[i],
          to = to,
          unit = "seconds"
        )
      # run only if channel exist
      channel <- wav[, channel.number]
      # get the frequency spectrum see seewave::spec for documentation
      spec <- myspec(channel, wav@samp.rate, fftw = T, plot = F)
      # get information like dominant frequency, etc.
      analysed <- specprop(spec)
      # get fundamental frequency
      #fnd <- fund(wav, fmax=2000, wl = 1024, plot = F)
      # only consider frequency range:
      filter <- spec[, 2][spec[, 1] < maxfreq / 1000]
      # take means:
      specmeans <-
        .colMeans(filter, length(filter) / maxfreq / binfreq, maxfreq / binfreq)
      frequencies <- seq(0, maxfreq - binfreq, binfreq)
      ampl.line <- paste(specmeans, collapse = ",")
      # write everything to csv
      write(
        paste(
          ampl.line,
          analysed$mean,
          analysed$sd,
          analysed$sem,
          analysed$median,
          analysed$mode,
          channel.number,
          #mean(fnd[,2], na.rm = T),
          as.character(current.timestamp),
          file,
          sep = ","
        ),
        file = csv,
        append = TRUE
      )
      current.timestamp <- current.timestamp + seconds
    }
    write(file, "passed.csv", append = T)
  }

main <-
  function(file,
           table,
           channel.number = 1,
           directory = "../../data/boden/") {
    print(file)
    library(tuneR)
    library(seewave)
    library(e1071)
    library(lubridate)
    tryCatch({
      if (substrRight(file, 3) == "wav") {
        analyse(file,
                directory,
                60 * 10,
                channel.number,
                paste(table, file, ".csv", sep = ""))
      }
    },
    error = function(cond) {
      write(paste(quoted(cond), quoted(file), sep = ","), "log.csv", append = T)
    },
    warning = function(cond) {
      write(paste(quoted(cond), quoted(file), sep = ","), "log.csv", append = T)
    },
    finally = {
    })
  }


directories <-
  c(
   "/media/diren/Bienenklänge BE A BEE/Daten_Vom Klang der Bienen/1.Speicherung13.4-18.4/miks/",
    "/media/diren/Bienenklänge BE A BEE/Daten_Vom Klang der Bienen/3.Speicherung bis 4.5/Miks/",
    "/media/diren/Bienenklänge BE A BEE/Daten_Vom Klang der Bienen/4.Speicherung 4.5-17.5/Miks/",
    "/media/diren/Bienenklänge BE A BEE/Daten_Vom Klang der Bienen/5.Speicherung 17.5-3.6/Miks/",
    "/media/diren/Bienenklänge BE A BEE/Daten_Vom Klang der Bienen/7.Speicherung 7.6-13.6/miks/",
    "/media/diren/Bienenklänge BE A BEE/Daten_Vom Klang der Bienen/8.Speicherung13.6-26.6/miks/"
  )
directories <- c("/home/diren/Documents/2018/beABee/data/mp3s/")

table <- "mp3s1/"
for (i in 1:length(directories)) {
  directory <- directories[i]
  files <- list.files(path = directory)
  print(files)
  
  # For a simple run:
  #files <- c("miks boden-03-may-20-52.wav", "miks boden-27-apr-12-46.mp3.wav")
  #lapply(files, main)
  # or:
  # main("miks boden-01-may-15-30.wav", "test.csv") # you might need to create an empty test.csv before
  
  # Run on one core only:
  #print(files)
  #lapply(files, main, table=table, directory=directory, channel.number=2)
  #
  # Run in parallel:
  #no_cores <- detectCores() # when I use to many cores, there's not enough memory, so I used max 3
  clust <- makeCluster(3)
  clusterExport(clust,
                list(
                  "substrRight",
                  "get.timestamp",
                  "analyse",
                  "myspec",
                  "quoted"
                ))
  parLapply(
    clust,
    files,
    main,
    table = table,
    directory = directory,
    channel.number = 2
  )
  stopCluster(clust)
}
# create Header-Line
fileConn <- file(paste(table, "aheader.csv", sep = ""))
frequencies <- seq(0, 1000 - 20, 20)
freq.line = paste(frequencies, sep = "", collapse = ", ")
writeLines(
  paste(
    freq.line,
    "mean",
    "sd",
    "sem",
    "median",
    "mode",
    "channel",
    "timestamp",
    "file",
    sep = ", "
  ),
  fileConn
)
close(fileConn)
#
# # Combine files after parallel run
# # use cat * > combined.csv in terminal, it is the fastes option!
#
# # Read the file, to test, if everything is there
# #combined <- read.table("tables/combined.csv", header=T, sep=",")
# # sort by timestamp
# #combined <- combined[order(combined$timestamp),]
# # plot fundamental frequency
# #plot(combined$fund)
