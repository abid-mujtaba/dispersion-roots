library(rjson)
library(ggplot2)

p <- ggplot()       # Initiate empty plot


# We have several series of data so we define a function that takes the file index (name/label) as well the linetype to use for that set and appends the appropriate plots
subplot <- function(p, index, ltype) {

    # Read data from json and csv files:

    VALUESFILE <- paste("value-", index, ".json", sep="")     # Create name of values file by concatenating with the PLOT label. sep="" removes the space between elements of the concatenation
    # To get the json data we read the file, concatenate all the lines and then use fromJSON() from the rjson library
    v <- fromJSON(paste(readLines(VALUESFILE), collapse="")) 

    DATAFILE <- paste("data-", index, ".csv", sep="")
    s <- read.csv(DATAFILE)


    # Create sbu-plots for the sequences
    for (seq in 1:7) {
        p <- p + geom_line(data=s[s$seq == seq,], aes(k_perp, omega), linetype=ltype)
    }


    # Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
    # and the use of paste() to combine normal text and symbols
    p <- p +
            labs(x = expression(k[perp]), y = expression(omega / omega[ce]))

    return (p)
}


# Repeatedly call subplot to add series of data
p <- subplot(p, "02-a", "solid")
p <- subplot(p, "02-b", "dashed")
p <- subplot(p, "02-c", "dotted")

#p <- p + xlim(0,10)        # Limit x-axis values


ggsave(file="plot-02.pdf", plot=p) 


# Source: http://www.cookbook-r.com/Graphs/Shapes_and_line_types/
