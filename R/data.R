########
## documentation for data in Vignette

#' Zachary's Karate Club graph
#'
#' Weighted adjacency matrix reporting interactions among
#' 34 nodes.
#'
#' @format a 34x34 matrix
#' @source package `igraphdata`
"adj_karate"

#' Zachary's Karate Club vertex faction assignment
#'
#' Vector reporting the assignment of nodes to communities.
#'
#' @format a 34-vector with the assignment of nodes to faction 1 or 2
#' @source package `igraphdata`
"vertexlabels"

#' Swiss MPs committee similarity matrix.
#'
#' **onlinesim_mat**: a similarity matrix of how similar two MPs are in their online
#' social media presence (shared supportees).
#'
#' @docType data
#'
#' @usage data(onlinesim_mat)
#'
#' @format 163x163 similarity matrix
#'
#' @keywords datasets
#'
# #' @references 
# #' (\href{}{})
#'
# #' @source \href{}{}
#'
"onlinesim_mat"

#' Swiss MPs committee affiliation data frame.
#'
#' **dtcommittee**: a list of committees each MP was part of during their stay in 
#' parliament
#'
#' @docType data
#'
#' @usage data(dtcommittee)
#'
#' @format 163x2 data.frame
#'
#' @keywords datasets
#'
# #' @references 
# #' (\href{}{})
#'
# #' @source \href{}{}
#'
"dtcommittee"

#' Swiss MPs attribute data frame.
#'
#' **dt**: contains different attributes of the 163 MPs, such as their names, 
#' their party affiliation (variable: *party*), their parliamentary group
#' affiliation (variable: *parlGroup*), the Canton (or state) they represent
#' (variable: *canton*), their gender  (variable: *gender*)
#' and date of birth  (variable: *birthdate*).
#'
#' @docType data
#'
#' @usage data(dt)
#'
#' @format 163x8 data.frame
#'
#' @keywords datasets
#'
# #' @references 
# #' (\href{}{})
#'
# #' @source \href{}{}
#'
"dt"

#' Swiss MPs network adjacency matrix
#'
#' **cospons_mat**: contains the adjacency matrix of 163 x 163 MPs.
#'
#' @docType data
#'
#' @usage data(cospons_mat)
#'
#' @format 163x163 adjacency matrix
#'
#' @keywords datasets
#'
# #' @references 
# #' (\href{}{})
#'
# #' @source \href{}{}
#'
"cospons_mat"

#' Highschool contact network multiplex representation
#'
#' **highschool.multiplex**: list containing the adjacency matrix of 327 x 327
#' highschool students, and the adjacency matrices correspodning to the 5
#' predictors used in Casiraghi2017.
#'
#' @docType data
#'
#' @usage data(highschool.multiplex)
#'
#' @format 6x327x327 list of adjacency matrices
#'
#' @keywords datasets
#'
#' @references Casiraghi, G. Multiplex Network Regression: How do relations
#'   drive interactions? 15 (2017).
#'
#'   Mastrandrea, R., Fournet, J. & Barrat, A. Contact patterns in a high
#'   school: A comparison between data collected using wearable sensors, contact
#'   diaries and friendship surveys. PLoS One 10, 1–26 (2015).
#'   
#' @source \href{sociopatterns.org}{https://sociopatters.org}
#'   
"highschool.multiplex"

#' Highschool contact network predictors
#'
#' **highschool.predictors**: list containing the adjacency matrices
#' corresponding to the 5 predictors used in Casiraghi2017.
#'
#' @docType data
#'
#' @usage data(highschool.predictors)
#'
#' @format 5x327x327 list of adjacency matrices
#'
#' @keywords datasets
#'
#' @references Casiraghi, G. Multiplex Network Regression: How do relations
#'   drive interactions? 15 (2017).
#'
#'   Mastrandrea, R., Fournet, J. & Barrat, A. Contact patterns in a high
#'   school: A comparison between data collected using wearable sensors, contact
#'   diaries and friendship surveys. PLoS One 10, 1–26 (2015).
#'
#' @source \href{sociopatterns.org}{https://sociopatters.org}
#'   
"highschool.predictors"

#' Highschool contact network adjacency matrix
#'
#' **contacts.adj**: contains the adjacency matrix of 327 x 327 highschool
#' students.
#'
#' @docType data
#'
#' @usage data(highschool.predictors)
#'
#' @format 327x327 adjacency matrix
#'
#' @keywords datasets
#'
#' @references Casiraghi, G. Multiplex Network Regression: How do relations
#' drive interactions? 15 (2017).
#'
#' Mastrandrea, R., Fournet, J. & Barrat, A. Contact patterns in a high school:
#' A comparison between data collected using wearable sensors, contact diaries
#' and friendship surveys. PLoS One 10, 1–26 (2015).
#' @source \href{https://sociopatterns.org}{https://sociopatters.org}
#'   
"contacts.adj"