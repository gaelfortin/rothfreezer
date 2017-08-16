# ---- View strain collection key ----
#' View of strain collection keyfile
#'
#' This view adds strain names and gene IDs from the strains table to the 
#' strain collection key.
#' 
#' @param src A connection to the rothfreezer database.
#' 
#' @return Returns a \code{tbl_sqlite} query, which can be further manipulated 
#' with \code{dplyr} verbs such as \link[dplyr]{filter}, and 
#' \link[dplyr]{mutate}. To execute your query and return a dataframe, use
#' \link[dplyr]{collect}. This view has the following fields:
#'
#' \item{\bold{strain_collection_id} (chr)}{ID of strain collection. 
#'   See the \code{strain_collection_info} table for a description of the 
#'   collection}
#' \item{\bold{strain_id} (chr)}{ID of strain. See the \code{strains} 
#'  table for more information about the strain}
#' \item{\bold{strain_name} (chr)}{Strain name (not necessarily unique).}
#' \item{\bold{gene_id} (chr)}{Gene ID (may cantain multiple gene IDs 
#'   separated by "|").}
#' \item{\bold{plate} (int)}{Plate number.}
#' \item{\bold{row} (int)}{Row numbr.}
#' \item{\bold{column} (int)}{Column number.}
#' \item{\bold{plate_control} (lgl)}{Is this strain a plate control in this 
#'   collection?}
#' 
#'
#' @examples \dontrun{
#' 
#' rad52 <- 
#'   view_strain_collection_keys() %>%
#'   filter(strain_name == 'rad52') %>%
#'   collect
#' }
#'
#' @export

view_strain_collection_keys <- function(src = src_rothfreezer()) {
  assert_that(is.src(src))
  left_join(
    tbl(src, 'strain_collections') %>% select_(~strain_collection_id:plate_control, ~excluded),
    tbl(src, 'strains') %>% select_(~strain_id, ~strain_name, ~gene_id),
    by = 'strain_id'
  ) %>%
  select_(~strain_collection_id, ~strain_id, ~strain_name, ~gene_id, ~plate:plate_control, ~excluded)
}
