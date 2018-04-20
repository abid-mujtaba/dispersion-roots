source('custom_themes.R')   # Contains custom_theme and custom_secondary_axis (for tick reflection)
source('subplots.R')

filename <- 'plot-02.tex'

p <- ggplot()

# Repeatedly call subplot to add series of data
p <- subplot(p, "02-a", "0.00")
p <- subplot(p, "02-b", "0.01")
p <- subplot(p, "02-c", "0.10")
p <- subplot(p, "02-d", "0.25")

# Add title, labels, and customize theme
#p <- p + ggtitle("Dispersion for $\\Lambda_h = 0$, $\\kappa_c = 3$, $\\kappa_h = 4$, $\\frac{n_{0h}}{n_{0e}} = 0.0$")
p <- p +
        scale_x_continuous("$k_\\perp \\rho_h$", sec.axis = custom_secondary_axis, expand  = c(0,0)) +
        scale_y_continuous("$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$", sec.axis = custom_secondary_axis) +
        scale_linetype_manual(name="$\\Lambda_c$", values=c("0.00"="solid", "0.01"="dotted", "0.10"="2282", "0.25"="dashed")) +
        expand_limits(x = 0) +    # Make graph start at x = 0 (snug up with left margin)
        custom_theme

if (plotTex) {

	write_plot(p, filename)

} else {

	write_plot(p, NULL, F)
}
