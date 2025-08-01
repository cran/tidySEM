#' National Identity, Discrimination and Depression
#'
#' These synthetic data are based on a study by Maene and colleagues,
#' which conducted an LCA with ordinal indicators on National, Regional,
#' and Heritage Identities in Flemish (Belgian) high-school students
#' with a migration background, and examined between class differences
#' in perceived discrimination by teachers and depressive symptoms.
#'
#' \tabular{lll}{
#'   \strong{Ethnic_1} \tab \code{ordered} \tab when I introduce myself,
#'   I would definitely say I belong to this group, answered on a 5-point
#'   Likert scale\cr
#'   \strong{Ethnic_2} \tab \code{ordered} \tab I have a strong sense of
#'   belonging to this group, answered on a 5-point Likert scale\cr
#'   \strong{Ethnic_3} \tab \code{ordered} \tab I see myself as a member of
#'   this group, answered on a 5-point Likert scale\cr
#'   \strong{Belgian} \tab \code{ordered} \tab Do you feel a member of
#'   the Belgian group, answered on a 10-point Likert scale\cr
#'   \strong{Flemish} \tab \code{ordered} \tab Do you feel a member of
#'   the Flemish group, answered on a 10-point Likert scale\cr
#'   \strong{age} \tab \code{numeric} \tab Participant age\cr
#'   \strong{sex} \tab \code{factor} \tab Participant sex\cr
#'   \strong{ses} \tab \code{numeric} \tab Socio-economic status,
#'   measured using the International Socio-Economic Index of Occupational
#'   Status (ISEI)\cr
#'   \strong{belgianborn} \tab \code{factor} \tab Whether or not the
#'   participant was born in Belgium\cr
#'   \strong{age_belgium} \tab \code{numeric} \tab Age at which the participant
#'   migrated to Belgium\cr
#'   \strong{vict_bully} \tab \code{factor} \tab Whether or not the
#'   participant has ever been the victim of peer bullying for any reason\cr
#'   \strong{vict_teacher} \tab \code{factor} \tab Whether or not the
#'   participant has ever been insulted, threatened, pushed, treated unfairly
#'   or excluded by teachers because of their foreign descent, language use,
#'   and skin colour\cr
#'   \strong{depression} \tab \code{numeric} \tab Scale scores of self-reported
#'   depressive feelings, assessed using the a ten-item scale with 5-point
#'   Likert-type response options
#' }
#' @docType data
#' @keywords datasets
#' @name maene_identity
#' @usage data(maene_identity)
#' @references Maene, C., D'hondt, F., Van Lissa, C. J., Thijs, J., & Stevens,
#' P. A. (2022). Perceived teacher discrimination and depressive feelings in
#' adolescents: the role of national, regional, and heritage identities in
#' Flemish schools. Journal of youth and adolescence, 51(12), 2281-2293.
#' \doi{10.1007/s10964-022-01665-7}
#' @format A data frame with 439 rows and 13 variables.
NULL
