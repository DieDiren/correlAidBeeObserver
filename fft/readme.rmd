---
title: "BeABee FFT"
subtitle: ""
author: "Diren"
---
The most import lines are the following, which are part of the analyse function.

```{r eval=F}
    # get the frequency spectrum see seewave::spec for documentation
    spec <- myspec(channel, wav@samp.rate, fftw = T, plot = F)
    # get information like dominant frequency, etc.
    analysed <- specprop(spec)
    # get fundamental frequency
    fnd <- fund(wav, fmax=2000, wl = 1024, plot = F)
    # only consider frequency range:
    filter <- spec[,2][spec[,1]<maxfreq/1000]
    # take means:
    specmeans <- .colMeans(filter, length(filter) / maxfreq/binfreq, maxfreq/binfreq)
```
myspec is the seewave:spec functions with some minor bugfix: [https://www.rdocumentation.org/packages/seewave/versions/1.0/topics/spec](https://www.rdocumentation.org/packages/seewave/versions/1.0/topics/spec)

specprob extracts the dominant frequency and stuff like that:
[https://www.rdocumentation.org/packages/seewave/versions/2.1.0/topics/specprop](https://www.rdocumentation.org/packages/seewave/versions/2.1.0/topics/specprop)

fund gives the fundamental frequency, but I am not sure, if this makes any sense:
[https://www.rdocumentation.org/packages/seewave/versions/2.1.0/topics/fund](https://www.rdocumentation.org/packages/seewave/versions/2.1.0/topics/fund)
