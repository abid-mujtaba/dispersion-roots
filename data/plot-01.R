source('custom_themes.R')

p <- ggplot()

# We have several series of data so we define a function that takes the file index (name/label) as well the linetype to use for that set and appends the appropriate plots
subplot <- function(p, index, parameter) {

    # Read data from csv file:

    DATAFILE <- paste("data-", index, ".csv", sep="")
    s <- read.csv(DATAFILE)


    # Create sub-plots for the sequences
    for (seq in 1:7) {
        ss <- s[s$seq == seq,]      # Create sub-set of the data for specified value of 'seq'
        p <- p + geom_line(data=ss, aes_(x=ss$k_perp, y=ss$omega, linetype=parameter))            # Use 'aes_' to gain access to local variable kappa_h (scope problems). This requires x= and y= to be declared explicitly. We set 'linetype' equal to the 'parameter' value and later manually provide a conversion from parameter value to the linetype
    }


    # Add title to plot. Note use of LaTeX expressions which are quoted and the back-slash is escaped
    # The .tex output will contain these LaTeX expressions inside the tikz diagram and they will be correctly rendered by pdflatex

    return (p)
}


# Repeatedly call subplot to add series of data
p <- subplot(p, "01-a", "0.00")
p <- subplot(p, "01-b", "0.01")
p <- subplot(p, "01-c", "0.10")
p <- subplot(p, "01-d", "0.25")

p <- p + ggtitle("Dispersion for $\\Lambda_c = 0$, $\\kappa_c = 2$, $\\kappa_h = 4$, $\\frac{n_{0h}}{n_{0e}} = 1.0$")
p <- p + scale_x_continuous("$k_\\perp$", sec.axis = custom_secondary_axis)
p <- p + scale_y_continuous("$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$", sec.axis = custom_secondary_axis)
p <- p + scale_linetype_manual(name="$\\Lambda_h$", values=c("0.00"="solid", "0.01"="dashed", "0.10"="dotted", "0.25"="dotdash"))
p <- p + custom_theme

write_plot(p, NULL, F)
