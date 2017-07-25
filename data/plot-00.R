library(rjson)
library(ggplot2)

PLOT <- "00-H"
VALUESFILE <- paste("value-", PLOT, ".json", sep="")     # Create name of values file by concatenating with the PLOT label. sep="" removes the space between elements of the concatenation
DATAFILE <- paste("data-", PLOT, ".csv", sep="")


# To get the json data we read the file, concatenate all the lines and then use fromJSON() from the rjson library
#v <- fromJSON(paste(readLines(VALUESFILE), collapse="")) 
s <- read.csv(DATAFILE)


p <- ggplot()       # Initiate empty plot

for (seq in 1:7) {
    p <- p + geom_line(data=s[s$seq == seq,], aes(k_perp, omega))
}


# Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
# and the use of paste() to combine normal text and symbols
p <- p +
         labs(x = expression(k[perp]), y = expression(omega / omega[ce]))

ggsave(file="plot-00.pdf", plot=p) 
