handle_request <- function(request) {
    command <- request$command %||% "health"
    result <- switch(command,
        health = health_result(),
        inspect_file = inspect_file(request),
        process_data = process_data(request),
        gsea = run_gsea(request),
        visualize = run_visualization(request),
        restore_state = restore_project_state(request),
        stop("Unknown analysis command: ", command)
    )
    emit(list(id = request$id, ok = TRUE, result = result))
}

run_worker_loop <- function(connection = file("stdin", open = "r", encoding = "UTF-8")) {
    on.exit(close(connection), add = TRUE)
    emit(list(event = "ready", status = "ready"))
    repeat {
        line <- readLines(connection, n = 1L, warn = FALSE)
        if (length(line) == 0L) break
        if (!nzchar(trimws(line[[1]]))) next
        request <- tryCatch(jsonlite::fromJSON(line[[1]], simplifyVector = FALSE), error = function(error) error)
        if (inherits(request, "error")) {
            emit_error(NULL, request)
            next
        }
        tryCatch(handle_request(request), error = function(error) emit_error(request$id %||% NULL, error))
    }
    invisible(TRUE)
}
