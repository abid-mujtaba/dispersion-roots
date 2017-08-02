library(rjson)
library(ggplot2)

p <- ggplot()       # Initiate empty plot


# We have several series of data so we define a function that takes the file index (name/label) as well the linetype to use for that set and appends the appropriate plots
subplot <- function(p, index, lambda) {

    # Read data from json and csv files:

    VALUESFILE <- paste("value-", index, ".json", sep="")     # Create name of values file by concatenating with the PLOT label. sep="" removes the space between elements of the concatenation
    # To get the json data we read the file, concatenate all the lines and then use fromJSON() from the rjson library
    v <- fromJSON(paste(readLines(VALUESFILE), collapse="")) 

    DATAFILE <- paste("data-", index, ".csv", sep="")
    s <- read.csv(DATAFILE)


    # Create sub-plots for the sequences
    for (seq in 1:7) {
        ss <- s[s$seq == seq,]      # Create sub-set of the data for specified value of 'seq'
        p <- p + geom_line(data=ss, aes_(x=ss$k_perp, y=ss$omega, linetype=lambda))            # Use 'aes_' to gain access to local variable kappa_h (scope problems). This requires x= and y= to be declared explicitly. We set 'linetype' equal to the 'kappa_h' value and later manually provide a conversion from kappa_h value to the linetype
    }


    # Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
    # and the use of paste() to combine normal text and symbols
    p <- p +
            labs(x = expression(k[perp]), y = expression(omega / omega[ce]))

    return (p)
}


# Repeatedly call subplot to add series of data
p <- subplot(p, "01-b-a", "0.00")
p <- subplot(p, "01-b-b", "0.10")
p <- subplot(p, "01-b-c", "0.20")
p <- p + scale_linetype_manual(name=expression(Lambda), values=c("0.00"="solid", "0.10"="dashed", "0.20"="dotted")) +      # The 'name' will be the title of the legend
         ggtitle(expression(paste("Roots of Dispersion Relation for ", kappa[c], " = 2, ", kappa[h], " = 4, ", frac(n[h0], n[e0]), " = 0.0")))

#p <- p + xlim(0,10)        # Limit x-axis values


ggsave(file="plot-01-b.pdf", plot=p) 


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
# Source for math expressions: http://vis.supstat.com/2013/04/mathematical-annotation-in-r/ 
