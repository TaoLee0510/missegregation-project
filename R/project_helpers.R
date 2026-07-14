raw_to_observations <- function(raw, dt) {
  steps <- names(raw)
  tags <- sort(unique(unlist(lapply(raw, names), use.names = FALSE)))
  x <- matrix(0L, nrow = length(tags), ncol = length(steps), dimnames = list(tags, steps))
  for (step in steps) x[names(raw[[step]]), step] <- as.integer(raw[[step]])
  list(x = x, dt = dt, passage_times = as.numeric(steps) * dt)
}

prepare_terminal_input <- function(observation_path) {
  observed <- readRDS(observation_path)
  final_path <- file.path(dirname(observation_path), "final_karyotypes.csv")
  final <- read.csv(final_path, stringsAsFactors = FALSE)
  if (!all(c("karyotype", "count") %in% names(final))) stop("Invalid final population for ", observation_path, call. = FALSE)
  terminal <- unique(final$karyotype[final$count > 0])
  if (!length(terminal)) stop("No nonzero endpoint karyotypes for ", observation_path, call. = FALSE)
  missing <- setdiff(terminal, rownames(observed$x))
  if (length(missing)) stop("Endpoint karyotypes are absent from the recorded trajectory: ", paste(missing, collapse = ", "), call. = FALSE)
  list(yi = list(x = observed$x[terminal, , drop = FALSE], dt = observed$dt), passage_times = observed$passage_times, terminal = terminal)
}

bind_status <- function(rows) {
  columns <- unique(unlist(lapply(rows, names), use.names = FALSE))
  do.call(rbind, lapply(rows, function(x) {
    x[setdiff(columns, names(x))] <- NA
    x[columns]
  }))
}
