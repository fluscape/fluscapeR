#' FluScape contact survey data (original)
#'
#' Household- and participant-level contact survey records from the original
#' FluScape study, with true (non-jittered) household coordinates.
#'
#' @format A data frame with 4710 rows and 14 variables:
#' \describe{
#'   \item{pid}{Participant identifier}
#'   \item{LOC_ID}{Location identifier}
#'   \item{HH_ID}{Household identifier}
#'   \item{PARTICIPANT_ID}{Participant identifier within household}
#'   \item{lat, long}{Latitude/longitude of the contact}
#'   \item{n.all, n.nothome, n.work, n.other}{Contact counts by category}
#'   \item{HH_Lat, HH_Long}{Latitude/longitude of the participant's household}
#'   \item{URBAN}{Indicator for urban location}
#'   \item{CHILD}{Indicator for whether the participant is a child}
#' }
"contacts_fluscape_V1"

#' FluScape contact survey data (jittered 100m)
#'
#' Same records as \code{\link{contacts_fluscape_V1}}, but with household
#' coordinates randomly jittered by up to 100m to allow public release
#' without disclosing exact household locations.
#'
#' @format A data frame with 4710 rows and 14 variables; see
#'   \code{\link{contacts_fluscape_V1}} for column definitions.
"contacts_V1_jittered_100m"
