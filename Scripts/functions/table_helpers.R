# Wrap long cell text before passing it to LaTeX tables.
wrap_table_text <- function(x, width = 58) {
  ifelse(is.na(x), NA_character_, stringr::str_wrap(as.character(x), width = width))
}

# Compact LaTeX table helper. Long tables use repeating headers, and wide tables
# are placed in landscape orientation when needed.
pdf_table <- function(x, caption, digits = 3, font_size = 8,
                      longtable = FALSE, scale_down = FALSE,
                      landscape = FALSE, widths = NULL) {
  names(x) <- gsub("_", " ", names(x), fixed = TRUE)
  table <- kableExtra::kbl(
    x,
    format = "latex",
    booktabs = TRUE,
    longtable = longtable,
    caption = caption,
    label = file_stub(caption),
    digits = digits,
    linesep = "",
    escape = TRUE
  )
  
  latex_options <- if (longtable) c("repeat_header") else c("HOLD_position")
  if (scale_down && !longtable) {
    latex_options <- c(latex_options, "scale_down")
  }
  
  table <- table |>
    kableExtra::kable_styling(
      latex_options = latex_options,
      font_size = font_size,
      full_width = FALSE,
      position = "center"
    )
  
  if (!is.null(widths)) {
    for (i in seq_along(widths)) {
      table <- table |>
        kableExtra::column_spec(i, width = widths[[i]])
    }
  }
  
  if (landscape) {
    table <- kableExtra::landscape(table)
  }
  table
}