set_up_dummy_vars <- function(res, comparisons){

  ###set up dummy variables
  for (i in 1:nrow(comparisons)) {
    if (!is.na(comparisons[i, "sample_type1"])) {
      res@colData[, as.character(comparisons[i, "name"])] <-
        ifelse(
          res@colData$cell_type == as.character(comparisons[i, "cell_type1"]) &
            res@colData$sample_name == as.character(comparisons[i, "sample_type1"]),
          paste(as.character(comparisons[i, "cell_type1"]),
                as.character(comparisons[i, "sample_type1"]), sep="_"),
          ifelse(
            res@colData$cell_type == as.character(comparisons[i, "cell_type2"]) &
              res@colData$sample_name == as.character(comparisons[i, "sample_type2"]),
            paste(as.character(comparisons[i, "cell_type2"]),
                  as.character(comparisons[i, "sample_type2"]), sep="_"),
            NA
          )
        )
    } else {
      res@colData[, as.character(comparisons[i, "name"])] <-
        ifelse(
          res@colData$cell_type == as.character(comparisons[i, "cell_type1"]),
          as.character(comparisons[i, "cell_type1"]),
          ifelse(
            res@colData$cell_type == as.character(comparisons[i, "cell_type2"]),
            as.character(comparisons[i, "cell_type2"]),
            NA
          )
        )
    }
  }
  return(res)
}
