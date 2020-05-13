## Changes from previous submission:
    
    # \dontrun{} should only be used if the example really cannot be executed
    # (e.g. because of missing additional software, missing API keys, ...) by the
    # user. That's why wrapping examples in \dontrun{} adds the comment ("# Not
    # run:") as a warning for the user. Does not seem necessary. Please unwrap the
    # examples if they are executable in < 5 sec, or create additionally small toy
    # examples to allow automatic testing. (You could also replace \dontrun{} with
    # \donttest, if it takes longer than 5 sec to be executed, but it would be
    # preferable to have automatic checks for functions. Otherwise, you can also
    # write some tests.)
    
  -> removed the existing \dontrun{} in favour of a \donttest{}. 
     Created faster examples when possible for the existing \donttest{}s and dropped them.

## Test environments
* local OS X install, R 3.6.3
* win-builder (devel, release and oldrelease)

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Giona Casiraghi <giona@ethz.ch>'

New submission

Possibly mis-spelled words in DESCRIPTION:
  Brandenberger (26:5)
  Casiraghi (19:33, 22:5, 23:5, 24:5, 25:5, 26:24, 27:5)
  ETH (21:84)
  Nanumyan (19:48, 22:20, 23:20, 25:20, 26:39)
  Scholtes (22:34, 23:34)
  arxiv (24:28)
  gHypEG (17:123)