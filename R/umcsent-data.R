#' Respondent-count data from Michigan Survey of Consumers
#'
#' Matrix of respondent counts. Rows are 479 months from January 1978 to November 2017. 
#' Columns are 1053 permutations of eight questions from the survey of consumers: 
#' PAGO, PEXP, PX1Q1, DUR,  BUS12, BUS5, UNEMP, and GOVT in the codebook. 
#' Each question has possible responses 1,3,5 and 8. This derived from the raw data 
#' by mapping response code 9 to 8, 2 to 1, and 4 to 5. Permutations with less than 38
#' respondents across all years are dropped to reduce 65,536 possible permutations to 
#' 1,053. 
#' 
#' @docType data
#'
#' @usage data(umcsent)
#'
#' @format A matrix with 479 rows and 1053 columns  
#' 
#' @keywords datasets
#'
#' @references Munro (2018) 
#' (\href{http://evanmunro.ca/files/discreteTS.pdf}{Working Paper})
#'
#' @source \href{https://data.sca.isr.umich.edu/sda-public/cgi-bin/hsda?harcsda+sca}{Survey of Consumers SDA Archive}
#'
#' @examples
#' data(umcsent)
"umcsent"