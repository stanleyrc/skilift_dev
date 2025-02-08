#' @name aggregate_events
#' @title aggregate_events
#' @description
#' Aggregate events per gene. 
#' Downstream of create_filtered_events().
#'
#' @param pair patient id
#' @return data.table or NULL
#' @export
aggregate_events <- function(
    pair,
    events_tbl,
	jabba_gg,
    cohort_type = "paired"
) {
	gg = readRDS(jabba_gg)
	karyotype_string = annotate_karyotype(gg)
}