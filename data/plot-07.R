library(rjson)
library(ggplot2)

p <- ggplot()       # Initiate empty plot


# We have several series of data so we define a function that takes the file index (name/label) as well the linetype to use for that set and appends the appropriate plots
subplot <- function(p, index, kappa_h) {

    # Read data from json and csv files:

    DATAFILE <- paste("data-", index, ".csv", sep="")
    s <- read.csv(DATAFILE)


    # Create sub-plots for the sequences
    for (seq in 1:7) {
        ss <- s[s$seq == seq,]      # Create sub-set of the data for specified value of 'seq'
        p <- p + geom_line(data=ss, aes_(x=ss$k_perp, y=ss$omega, linetype=kappa_h))            # Use 'aes_' to gain access to local variable kappa_h (scope problems). This requires x= and y= to be declared explicitly. We set 'linetype' equal to the 'kappa_h' value and later manually provide a conversion from kappa_h value to the linetype
    }


    # Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
    # and the use of paste() to combine normal text and symbols
    p <- p +
            labs(x = expression(k[perp]), y = expression(omega / omega[ce]))

    return (p)
}


# Repeatedly call subplot to add series of data
p <- subplot(p, "07-a", "1.6")
p <- subplot(p, "07-b", "4.0")
p <- subplot(p, "07-c", "10.0")
p <- p + scale_linetype_manual(name=expression(kappa[h]), values=c("1.6"="solid", "4.0"="dashed", "10.0"="dotted")) +      # The 'name' will be the title of the legend
         ggtitle(expression(paste("Roots of Dispersion Relation for ", Lambda, " = 0.00, ", kappa[c], " = 2, ", frac(n[h0], n[e0]), " = 1.0")))

#p <- p + xlim(0,10)        # Limit x-axis values


ggsave(file="plot-07.pdf", plot=p) 


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
# Source for math expressions: http://vis.supstat.com/2013/04/mathematical-annotation-in-r/ 