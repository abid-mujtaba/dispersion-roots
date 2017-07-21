library(ggplot2)

s <- read.csv("data.csv")


p <- ggplot()       # Initiate empty plot

for (seq in 1:7) {
    p <- p + geom_line(data=s[s$seq == seq,], aes(k_perp, omega))
}

# Add title to plot. Note use of expression() to access the plotmath ability to typecast equations/sybmols
# and the use of paste() to combine normal text and symbols
p <- p +
         ggtitle(expression(paste("Roots of VC Dispersion Relation (", Lambda, " = 0.2)"))) +
         labs(x = expression(k[perp]), y = expression(omega / omega[ce]))

ggsave(file="plot.pdf", plot=p) 
