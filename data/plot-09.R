source('custom_themes.R')   # Contains custom_theme and custom_secondary_axis (for tick reflection)
source('subplots.R')

filename <- 'plot-09.tex'

p <- ggplot()


# Repeatedly call subplot to add series of data
p <- subplot(p, "09-a", "10")
p <- subplot(p, "09-b", "100")
p <- subplot(p, "09-c", "1000")

# Add title, labels, and customize theme
p <- p + ggtitle("Dispersion for $\\Lambda_c = 0.01$, $\\Lambda_h = 0.1$, $\\kappa_c = 2$, $\\kappa_h = 4$, $\\displaystyle \\frac{n_{0h}}{n_{0e}} = 0.5$")
p <- p +
        scale_x_continuous("$k_\\perp$", sec.axis = custom_secondary_axis, expand  = c(0,0)) +
        scale_y_continuous("$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$", sec.axis = custom_secondary_axis) +
        scale_linetype_manual(name="$\\displaystyle \\frac{T_h}{T_c}$", values=c("10"="solid", "100"="dashed", "1000"="dotted")) +
        expand_limits(x = 0) +    # Make graph start at x = 0 (snug up with left margin)
        custom_theme

if (plotTex) {

	write_plot(p, filename)

} else {

	write_plot(p, NULL, F)
}
