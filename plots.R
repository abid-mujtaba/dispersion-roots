library(ggplot2)

s <- read.csv("data.csv")

s1 <- s[s$seq == 1,]
plot(s1$k_perp, s1$omega, "l", xlim=c(0,100), ylim=c(1,8))

for (seq in 2:7) {
        ss <- s[s$seq == seq,]
        lines(ss$k_perp, ss$omega, "l")
}
