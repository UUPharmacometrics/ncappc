# Remove CRAN note on no visible binding for global variable
utils::globalVariables(c('.'))

# ignore unused imports
ignore_unused_imports <- function() {
  bookdown::html_document2()
}
