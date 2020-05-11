## Changes from previous submission:
    
    # Please shorten the title to a maximum of 65 characters.
    # Acronyms/Abbreviations can be used on their own in the title 
    # as long as they are explained in the description field.
  
  -> Title shortened to 63 characters.

    # Please do not ship the full license but only additional restrictions,
    # otherwise omit the file and its reference entirely. 
    # A copy of the GPL is part of R.
  
  -> License file removed together with its reference in the DESCRIPTION.

    # Missing Rd-tags:
    #      homophily_stat.Rd: \value
    #      print.bccm.Rd: \value
    #      print.ghype.Rd: \value
    #      print.nrm.Rd: \value
    #      print.nrm.selection.Rd: \value
    #      summary.nrm.Rd: \value
    #      summary.nrm.selection.Rd: \value
  
  -> added \value or specified a @describeIn for all relevant functions.

    # You write information messages to the console that 
    # cannot be easily suppressed.
    # It is more R like to generate objects that can be used 
    # to extract the information a user is interested in, and 
    # then print() that object.
    # Instead of print()/cat() rather use message()/warning() 
    # or if(verbose)cat(..) if you really have to write text to the console.
    # (except for print() and summary() functions)
    
  -> removed all print()/cat() in favour of message().

    # Please always add all authors, contibutors and copyright holders 
    # in the Authors@R field with the appropriate roles.
    # From the CRAN policies you agreed to:
    # "The ownership of copyright and intellectual property rights of all
    # components of the package must be clear and unambiguous (including 
    # from the authors specification in the DESCRIPTION file). 
    # Where code is copied (or derived) from the work of others 
    # (including from R itself), care must be taken that any
    # copyright/license statements are preserved and authorship 
    # is not misrepresented.
    # Preferably, an ‘Authors@R’ would be used with ‘ctb’ roles for the
    # authors of such code. Alternatively, the ‘Author’ field should
    # list these authors as contributors.
    # Where copyrights are held by an entity other than the package authors,
    # this should preferably be indicated via ‘cph’ roles in the ‘Authors@R’
    # field, or using a ‘Copyright’ field (if necessary referring to an
    # inst/COPYRIGHTS file)."
    #     Missing: GV
  
  -> Added missing reference to GV in Authors@R.

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