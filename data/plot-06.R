source('custom_themes.R')   # Contains custom_theme and custom_secondary_axis (for tick reflection)
source('subplots.R')

filename <- 'plot-06.tex'

p <- ggplot()


# Repeatedly call subplot to add series of data
p <- subplot(p, "06-a", "2.6")
p <- subplot(p, "06-b", "3.0")
p <- subplot(p, "06-c", "$\\inf$")


# Add title, labels, and customize theme
#p <- p + ggtitle("Dispersion for $\\Lambda_c = 0.01$, $\\Lambda_h = 0.1$, $\\kappa_h = 4$, $\\displaystyle \\frac{n_{0h}}{n_{0e}} = 0.5$, $\\displaystyle \\frac{T_h}{T_c} = 101.695$")
p <- p +
        scale_x_continuous("$k_\\perp \\rho_h$", sec.axis = custom_secondary_axis, expand  = c(0,0)) +
        scale_y_continuous("$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$", sec.axis = custom_secondary_axis) +
        scale_linetype_manual(name="$\\kappa_c$", breaks=c("2.6", "3.0", "$\\inf$"), values=c("solid", "dashed", "dotted")) +
        expand_limits(x = 0) +    # Make graph start at x = 0 (snug up with left margin)
        custom_theme

write_plot(p, filename, texOutput = T)
