# Shared typography contract for every ATOP figure.

atop_axis_linewidth <- 0.25

atop_axis_line <- function() {
    ggplot2::element_line(colour = "black", linewidth = atop_axis_linewidth, lineend = "butt")
}

atop_font_size_mm <- function(font_size_pt) {
    as.numeric(font_size_pt) / (72.27 / 25.4)
}

atop_plot_text_theme <- function(font_size_pt, font_family = "Helvetica") {
    ggplot2::theme(
        axis.line = atop_axis_line(),
        axis.ticks = atop_axis_line(),
        text = ggplot2::element_text(
            family = font_family,
            size = font_size_pt,
            colour = "#000000"
        ),
        plot.title = ggplot2::element_text(
            family = font_family,
            size = font_size_pt,
            hjust = 0.5,
            face = "plain",
            colour = "#000000",
            margin = ggplot2::margin(0, 0, 4, 0)
        ),
        axis.text = ggplot2::element_text(size = font_size_pt, colour = "#000000"),
        axis.title = ggplot2::element_text(size = font_size_pt, colour = "#000000"),
        legend.title = ggplot2::element_text(size = font_size_pt, colour = "#000000"),
        legend.text = ggplot2::element_text(size = font_size_pt, colour = "#000000"),
        strip.text = ggplot2::element_text(size = font_size_pt, colour = "#000000")
    )
}
