
























create_global_permutation_tracker <- function(total_permutations, progress_callback = NULL) {
    start_time <- Sys.time()
    completed_permutations <- 0
    last_update_time <- start_time

    update_global_progress <- function(perm_increment = 0, custom_message = NULL, force_update = FALSE) {
        current_time <- Sys.time()
        completed_permutations <<- completed_permutations + perm_increment
        elapsed_time <- as.numeric(current_time - start_time, units = "secs")


        time_since_last_update <- as.numeric(current_time - last_update_time, units = "secs")
        if (!force_update && time_since_last_update < 0.5 && perm_increment > 0) {
            return(NULL)
        }
        last_update_time <<- current_time

        if (completed_permutations > 0 && elapsed_time > 0) {
            progress_ratio <- min(completed_permutations / total_permutations, 1.0)


            perms_per_sec <- completed_permutations / elapsed_time
            remaining_perms <- max(0, total_permutations - completed_permutations)
            estimated_remaining_time <- if (remaining_perms > 0 && perms_per_sec > 0) remaining_perms / perms_per_sec else 0


            format_time <- function(seconds) {
                if (seconds < 60) {
                    paste0(round(seconds, 0), "s")
                } else if (seconds < 3600) {
                    paste0(round(seconds / 60, 1), "m")
                } else {
                    paste0(round(seconds / 3600, 1), "h")
                }
            }


            format_speed <- function(speed) {
                if (speed >= 1000) {
                    paste0(round(speed / 1000, 1), "K perms/sec")
                } else if (speed >= 1) {
                    paste0(round(speed, 1), " perms/sec")
                } else if (speed > 0) {
                    paste0(round(60 / speed, 1), " sec/perm")
                } else {
                    "Calculating..."
                }
            }


            if (is.null(custom_message)) {
                speed_text <- format_speed(perms_per_sec)
                progress_text <- paste0("Permutation Test: ", completed_permutations, "/", total_permutations)

                if (estimated_remaining_time > 1 && remaining_perms > 0) {
                    eta_text <- format_time(estimated_remaining_time)
                    message <- paste0(progress_text, " (", speed_text, ") ⏱️ ETA: ", eta_text)
                } else {
                    message <- paste0(progress_text, " (", speed_text, ")")
                }
            } else {
                message <- custom_message
            }


            formatted_message <- paste0(message, " [", round(progress_ratio * 100, 1), "%]")
            if (is.function(progress_callback)) {
                progress_callback(formatted_message)
            } else {
                cat(paste0("\r", formatted_message))
                if (completed_permutations >= total_permutations) cat("\n")
            }

            return(list(
                completed = completed_permutations,
                total = total_permutations,
                progress = progress_ratio,
                elapsed_time = elapsed_time,
                estimated_remaining = estimated_remaining_time,
                speed = perms_per_sec
            ))
        }

        return(NULL)
    }

    return(update_global_progress)
}
