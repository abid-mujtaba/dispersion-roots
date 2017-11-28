source('custom_themes.R')   # Contains custom_theme and custom_secondary_axis (for tick reflection)
source('subplots.R')

p <- ggplot()


# Repeatedly call subplot to add series of data
p <- subplot(p, "07-a", "0.0")
p <- subplot(p, "07-b", "0.3")
p <- subplot(p, "07-c", "0.5")
p <- subplot(p, "07-d", "0.8")
p <- subplot(p, "07-e", "1.0")


# Add title, labels, and customize theme
p <- p + ggtitle("Dispersion for $\\Lambda_c = 0.01$, $\\Lambda_h = 0.1$, $\\kappa_c = 2$, $\\kappa_h = 4$, $\\displaystyle \\frac{T_h}{T_c} = 101.695$")
p <- p +
        scale_x_continuous("$k_\\perp$", sec.axis = custom_secondary_axis, expand  = c(0,0)) +
        scale_y_continuous("$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$", sec.axis = custom_secondary_axis) +
        scale_linetype_manual(name="$\\displaystyle \\frac{n_{0h}}{n_{0e}}$", values=c("0.0"="solid", "0.3"="dashed", "0.5"="dotted", "0.8"="dotdash", "1.0"="longdash")) +
        expand_limits(x = 0) +    # Make graph start at x = 0 (snug up with left margin)
        custom_theme

# write_plot(p, 'plot-07.tex')
write_plot(p, NULL, F)
