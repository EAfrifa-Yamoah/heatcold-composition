# 00_harvest.R
# Copy verbatim GBD 2023 connector results (saved by the session harness as text files) into
# data/raw/ under canonical names derived from the content of each file (risk, cause, measure,
# age basis, year), so that no file name depends on the order in which calls were made.
# Usage: Rscript R/00_harvest.R <tool-results directory>
args <- commandArgs(trailingOnly = TRUE)
src <- if (length(args)) args[1] else "/root/.claude/projects/-home-claude/4a32cf05-c4a9-5296-84b9-e2caa28798ef/tool-results"
files <- list.files(src, pattern = "^mcp-GBD_2023_Collab-.*\\.txt$", full.names = TRUE)
slug <- function(x) gsub("_+$", "", gsub("[^a-z0-9]+", "_", tolower(x)))
log <- list()
for (f in files) {
  txt <- readLines(f, warn = FALSE)
  fmt <- "csv"; x <- NULL
  if (length(txt) && substr(txt[1], 1, 1) == "{") {
    j <- jsonlite::fromJSON(paste(txt, collapse = "\n"))
    if (!is.null(j$rows)) {            # JSON response: values rounded to GBD publication precision
      if (!is.null(j$truncation)) next # trimmed responses are never used
      x <- j$rows; fmt <- "json"; tc <- textConnection("txt2", "w", local = TRUE); write.csv(x, tc, row.names = FALSE); close(tc); txt <- txt2
    } else txt <- strsplit(j$result, "\n")[[1]]
  }
  if (is.null(x)) x <- tryCatch(read.csv(text = paste(txt, collapse = "\n"), check.names = FALSE), error = function(e) NULL)
  if (is.null(x) || !nrow(x)) next
  if ("note" %in% names(x) && any(grepl("TRUNCATED", x$note))) next
  tool <- sub("^mcp-GBD_2023_Collab-(.*)-[0-9]+\\.txt$", "\\1", basename(f))
  nm <- names(x); u <- function(col) if (col %in% nm) paste(unique(x[[col]]), collapse = "+") else NA
  yr <- if ("year" %in% nm) { y <- sort(unique(x$year)); if (length(y) == 1) y else paste0(min(y), "to", max(y)) } else NA
  ag <- if ("age_group" %in% nm && length(unique(x$age_group)) > 3) "byage" else slug(u("age_group"))
  ag <- sub("age_standardized", "as", sub("all_ages", "all", ag))
  key <- if (tool == "get_covariate") sprintf("cov__%s__%s__n%d", slug(u("covariate_name")), yr, length(unique(x$location_id)))
    else if ("rei_id" %in% nm) sprintf("%s__r%s__c%s__%s__%s__%s__n%d", sub("get_", "", tool), u("rei_id"), if ("cause_id" %in% nm) u("cause_id") else slug(u("cause_name")), slug(u("measure")), ag, yr, length(unique(x$location_id)))
    else if ("risk_name" %in% nm) sprintf("%s__%s__%s__%s__%s__%s__n%d", sub("get_", "", tool), slug(u("risk_name")), slug(u("cause_name")), slug(u("measure")), ag, yr, length(unique(x$location_id)))
    else sprintf("%s__c%s__%s__%s__%s__n%d", sub("get_", "", tool), if ("cause_id" %in% nm) u("cause_id") else slug(u("cause_name")), slug(u("measure")), ag, yr, length(unique(x$location_id)))
  if (tool == "get_covariate") key <- paste0(key, "__", substr(slug(x$location_name[1]), 1, 12))
  if (tool == "get_multi_location_trend" && !(1 %in% x$location_id)) key <- paste0(key, "__", substr(slug(x$location_name[1]), 1, 12))
  if (fmt == "json") key <- paste0(key, "__jsonrounded")
  writeLines(txt, file.path("data/raw", paste0(key, ".csv")))
  log[[length(log) + 1]] <- data.frame(file = paste0(key, ".csv"), rows = nrow(x), format = fmt, source = basename(f))
}
log <- do.call(rbind, log); log <- log[!duplicated(log$file, fromLast = TRUE), ]
write.csv(log[order(log$file), ], "data/raw/_harvest_log.csv", row.names = FALSE)
cat("harvested", nrow(log), "distinct extracts\n")
