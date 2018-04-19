# R script containing utility functions for plotting
# These functions are used for inverting and mirror ticks

# Source: http://stackoverflow.com/a/29023682/2926226 (with minor edits)

library(tikzDevice)
library(ggplot2)

# Function for inverting the ticks (point inwards) on a graph
# Personal edit: Removed first arg to theme()
invert_ticks <- function() {

        t <- theme(
                     axis.ticks.length=unit(-0.2, "cm"),
                     axis.text.x = element_text(margin = margin(t = .5, unit = "cm")),
                     axis.text.y = element_text(margin = margin(r = .5, unit = "cm"))
                  )

        return (t)
}


# Define the thematic elements common to each plot. Use the black-and-white theme with inverted ticks and increased text sizes
custom_theme  <- theme_bw() +
        invert_ticks() +
        theme(
            plot.title = element_text(hjust = 0.5, size=16, margin = margin(b = 30, unit="pt")),    # Increase size of title and shift it up
            axis.title.y = element_text(size=12, angle=0, vjust=0.5),              # Make the label upright
            axis.title.x = element_text(size=12),
            axis.text = element_text(size=10, color="black"),
            legend.title = element_text(size=12),
            legend.text = element_text(size=10),
            axis.line = element_line(size=0.5),     # Increase thickness of axis line and ticks
            axis.ticks = element_line(size=0.5),
            panel.grid.major = element_blank(),     # Turn off all grid lines. Reactive y-grid explicitly below.
            panel.grid.minor = element_blank(),
            panel.grid.major.y = element_line(color = "grey70", size=0.4),    # Make major grid more visible.
            panel.grid.minor.y = element_line(color = "grey70", size=0.4)) +                           # Turn off minor grid.
        theme(plot.margin = margin(t = 10, b = 10, l = 10, r = 20, unit="points"))              # Needed to fix margins messed with by scale_y_continous (sy)

# Define axis duplication for secondary axis
custom_secondary_axis <- sec_axis(~ ., labels=NULL)


write_plot <- function(plot, filename, use_tikz = T) {
        if (use_tikz) {
                # Set tikz output device which will be used by the next 'print' command.
                tikz(file=filename, width = 7, height = 7)
                print(plot)
                dev.off()               # Without this the file won't be written to and closed properly
        } else {
                print(plot)
        }
}
