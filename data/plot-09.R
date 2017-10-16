library(rjson)
library(tikzDevice)
library(ggplot2)


# Open the tikz device for output
tikz(file='plot-09.tex', width=7, height=7)

p <- ggplot() + theme_bw()       # Initiate empty plot


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
p <- subplot(p, "09-a", "10")
p <- subplot(p, "09-b", "100")
p <- subplot(p, "09-c", "1000")

p <- p + labs(x = "$k_\\perp$", y = "$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$") +
         ggtitle("Dispersion for $\\Lambda_c = 0.01$, $\\Lambda_h = 0.1$, $\\kappa_c = 2$, $\\kappa_h = 4$, $\\displaystyle \\frac{n_{0h}}{n_{0e}} = 0.5$")

p <- p + scale_linetype_manual(name="$\\displaystyle \\frac{T_h}{T_c}$", values=c("10"="solid", "100"="dashed", "1000"="dotted"))

# Rotate the y-axis title and make its vertical justification centered
p <- p + theme(axis.title.y = element_text(angle = 0, vjust=0.5))

axis.title = element_text(size = rel(1.5))      # Increase size of axis titles


print(p)        # This will create the output using tikzDevice
dev.off()       # Necessary to close tikzDevice so the .tex file is written


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
# Source for tikzdevice: http://iltabiai.github.io/tips/latex/2015/09/15/latex-tikzdevice-r.html
