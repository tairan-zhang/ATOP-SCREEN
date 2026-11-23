# Shared conventions for gene-level figures. Plot builders own their data
# transformation; this module only provides visual and validation primitives.

atop_require_columns <- function(data, required, context) {
    missing <- setdiff(required, names(data))
    if (length(missing) > 0L) {
        stop(context, " is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
    }
}
