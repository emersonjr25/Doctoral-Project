##### FUNCTION:CREATING FOLDER TO REPLICATION #####
save_traits <- function() {
  save_extract("traits", rep)
}

save_extract <- function(element, replicate_number) {
  #browser()
  config <- dynGet("config")
  data <- dynGet("data")
  vars <-  dynGet("vars")
  save_landscape()
  folder_per_rep <- paste0('traits', replicate_number)
  dir.create(file.path(config$directories$output, folder_per_rep), showWarnings=FALSE, recursive=TRUE)
  tmp <- lapply(data$all_species, function(x){return(x[[element]])})
  names(tmp) <- sapply(data$all_species, function(x){x$id})
  saveRDS(object = tmp,
          file = file.path(config$directories$output, folder_per_rep, paste0(element, "_t_", vars$ti, ".rds")))
}
