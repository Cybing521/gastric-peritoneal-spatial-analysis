args <- commandArgs(trailingOnly = TRUE)
rda <- args[[1]]
out <- args[[2]]
nm <- load(rda)
obj <- get(nm[1])
inter <- NULL
if (is.data.frame(obj) || is.matrix(obj)) {
  inter <- as.data.frame(obj)
} else if (is.list(obj)) {
  if (!is.null(obj$interaction)) inter <- as.data.frame(obj$interaction)
}
if (is.null(inter)) stop("No interaction table in ", paste(nm, collapse = ","))
cols <- colnames(inter)
message("columns: ", paste(cols, collapse = " | "))
need <- c("ligand", "receptor")
if (!all(need %in% cols)) stop("missing ligand/receptor")
keep <- inter[, intersect(c("ligand", "receptor", "pathway_name", "interaction_name"), cols), drop = FALSE]
write.csv(keep, out, row.names = FALSE)
message("wrote ", nrow(keep), " rows")
