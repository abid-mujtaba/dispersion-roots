library(rjson)
library(ggplot2)

p <- ggplot()       # Initiate empty plot


# We have several series of data so we define a function that takes the file index (name/label) as well the linetype to use for that set and appends the appropriate plots
subplot <- function(p, index, parameter) {

    # Read data from csv files:
    DATAFILE <- paste("data-", index, ".csv", sep="")
    s <- read.csv(DATAFILE)


    # Create sub-plots for the sequences
    for (seq in 1:7) {
        ss <- s[s$seq == seq,]      # Create sub-set of the data for specified value of 'seq'
        p <- p + geom_line(data=ss, aes_(x=ss$k_perp, y=ss$omega, linetype=parameter))            # Use 'aes_' to gain access to local variable 'parameter' (scope problems). This requires x= and y= to be declared explicitly. We set 'linetype' equal to the 'parameter' value and later manually provide a conversion from parameter value to the linetype
    }


    # Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
    # and the use of paste() to combine normal text and symbols
    p <- p +
            labs(x = expression(k[perp]), y = expression(omega / omega[ce]))

    return (p)
}


# Repeatedly call subplot to add series of data
p <- subplot(p, "09-a", "0.0")
p <- subplot(p, "09-b", "0.3")
p <- subplot(p, "09-c", "0.5")
p <- subplot(p, "09-d", "0.8")
p <- subplot(p, "09-e", "1.0")
p <- p + scale_linetype_manual(name=expression(frac(n[h0], n[e0])), values=c("0.0"="solid", "0.3"="dotted", "0.5"="dashed", "0.8"="dotdash", "1.0"="longdash")) +      # The 'name' will be the title of the legend
         ggtitle(expression(paste("Roots of Dispersion Relation for ", kappa[c], " = 2, ", kappa[h], " = 4, ", Lambda, " = 0.2")))

#p <- p + xlim(0,10)        # Limit x-axis values


ggsave(file="plot-09.pdf", plot=p)


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
# Source for math expressions: http://vis.supstat.com/2013/04/mathematical-annotation-in-r/