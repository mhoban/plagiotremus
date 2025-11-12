library(R6)
library(worrms)

Taxa <- R6Class(
  "Taxa",
  private = list(
    taxa = NULL
  ),
  public = list(
    initialize = function(genus = '') {
      id <- wm_name2id(genus)
      private$taxa <- wm_record(id) %>%
        bind_rows() %>%
        bind_rows(
          wm_children(id)
        )
    },
    cite = function(sp,abbrev=FALSE,italic=TRUE) {
      it <- ifelse(italic,"*","")
      tt <- private$taxa %>%
        filter(scientificname == sp)
      if (nrow(tt) > 0) {
        tt %>% 
          mutate(
            newname = if_else(
              abbrev,
              str_replace(scientificname,"([A-Z])\\w+ (\\w+)","\\1. \\2"),
              scientificname
            ),
            citation = str_glue("{it}{newname}{it} {authority}")
          ) %>%
          pull(citation) 
      } else {
        str_glue("{it}{sp}{it}")
      }
    },
    get_taxa = function() private$taxa
  )
)