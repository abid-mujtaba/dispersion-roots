library(rjson)
library(tikzDevice)
library(ggplot2)


# Open the tikz device for output
tikz(file='plot-07.tex', width=7, height=7)

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
p <- subplot(p, "07-a", "0.0")
p <- subplot(p, "07-b", "0.3")
p <- subplot(p, "07-c", "0.5")
p <- subplot(p, "07-d", "0.8")
p <- subplot(p, "07-e", "1.0")

p <- p + labs(x = "$k_\\perp$", y = "$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$") +
         ggtitle("Dispersion for $\\Lambda_c = 0.01$, $\\Lambda_h = 0.1$, $\\kappa_c = 2$, $\\kappa_h = 4$, $\\displaystyle \\frac{T_h}{T_c} = 101.695$")

p <- p + scale_linetype_manual(name="$\\displaystyle \\frac{n_{0h}}{n_{0e}}$", values=c("0.0"="solid", "0.3"="dashed", "0.5"="dotted", "0.8"="dotdash", "1.0"="longdash"))

# Alter elements of the plot
p <- p + theme(axis.title.y = element_text(angle = 0, vjust=0.5),       # Rotate the y-axis title and make its vertical justification centered
               plot.title = element_text(margin = margin(t=10, b=15, l=0, r=0, unit="pt")),           # Increase vertical margin of title to adjust for fractions
               legend.title = element_text(margin = margin(t=10, b=35, l=0, r=0, unit="pt")))         # Increase vertical margin of legend title to adjust for fractions (Does NOT work - BUG)


print(p)        # This will create the output using tikzDevice
dev.off()       # Necessary to close tikzDevice so the .tex file is written


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
# Source for tikzdevice: http://iltabiai.github.io/tips/latex/2015/09/15/latex-tikzdevice-r.html
