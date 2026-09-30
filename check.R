## Check ctv is correct

library(ctv)

ctv <- "TimeSeries"
ctvfile <- here::here(paste0(ctv, ".md"))
htmlfile <- here::here(paste0(ctv, ".html"))

# Set CRAN mirror
r <- getOption("repos")
r["CRAN"] <- "https://cloud.r-project.org"
options(repos = r)

# Run the check
check_ctv_packages(ctvfile)

# Create html file from ctv file
ctv2html(read.ctv(ctvfile), htmlfile)

# Check the contents list: every link must match a heading anchor,
# and every heading (other than Links) must appear in the contents
html <- paste(readLines(htmlfile, warn = FALSE), collapse = "\n")
ids <- regmatches(html, gregexpr('<h[34] id="[^"]+"', html))[[1]]
ids <- sub('.*id="([^"]+)"', "\\1", ids)
md <- readLines(ctvfile, warn = FALSE)
toc <- unlist(regmatches(md, gregexpr("\\]\\(#[^)]+\\)", md)))
toc <- gsub("^\\]\\(#|\\)$", "", toc)
broken <- setdiff(toc, ids)
unlisted <- setdiff(ids, c(toc, "contents"))
if (length(broken)) {
  cat("Broken contents links:\n", paste0("  #", broken, "\n"), sep = "")
}
if (length(unlisted)) {
  cat(
    "Headings missing from contents:\n",
    paste0("  #", unlisted, "\n"),
    sep = ""
  )
}
if (!length(broken) && !length(unlisted)) {
  cat("Contents links OK.\n")
}

# View html file
browseURL(htmlfile)

cat("Done.\n")
