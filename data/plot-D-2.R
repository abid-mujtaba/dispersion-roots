library(rjson)
library(ggplot2)

p <- ggplot()       # Initiate empty plot


# We have several series of data so we define a function that takes the file index (name/label) as well the linetype to use for that set and appends the appropriate plots
subplot <- function(p, index, parameter) {

    # Read data from csv files:
    DATAFILE <- paste("data-", index, ".csv", sep="")
    s <- read.csv(DATAFILE)


    p <- p + geom_line(data=s, aes_(x=s$omega, y=s$D, linetype=parameter))            # Use 'aes_' to gain access to local variable 'parameter' (scope problems). This requires x= and y= to be declared explicitly. We set 'linetype' equal to the 'parameter' value and later manually provide a conversion from parameter value to the linetype


    # Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
    # and the use of paste() to combine normal text and symbols
    p <- p +
            labs(x = expression(omega), y = expression(D))
        

    return (p)
}


# Repeatedly call subplot to add series of data
p <- subplot(p, "D-2-a", "1.6")
p <- subplot(p, "D-2-b", "2.0")
p <- subplot(p, "D-2-c", "2.4")
p <- p + scale_linetype_manual(name=expression(kappa[c]), values=c("1.6"="solid", "2.0"="dotted", "2.4"="dashed")) +      # The 'name' will be the title of the legend
         ggtitle(expression(paste("Roots of Dispersion Relation for ", frac(n[h0], n[e0]), " = 0.0, ", Lambda, " = 0.00 ", k[perp], " = 5.0")))

#p <- p + xlim(0,10)        # Limit x-axis values
p <- p + ylim(-50,50)


ggsave(file="plot-D-2.pdf", plot=p) 


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
# Source for math expressions: http://vis.supstat.com/2013/04/mathematical-annotation-in-r/ 
