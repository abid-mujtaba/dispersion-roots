source('custom_themes.R')   # Contains custom_theme and custom_secondary_axis (for tick reflection)
source('subplots.R')

filename <- 'plot-12.tex'

p <- ggplot()


# Repeatedly call subplot to add series of data
p <- subplot(p, "12-a", "Maxwellian")
p <- subplot(p, "12-b", "Double-Kappa")
p <- subplot(p, "12-c", "Cairns")
p <- subplot(p, "12-d", "VC")

# Add title, labels, and customize theme
#p <- p + ggtitle("Dispersion for $\\displaystyle \\frac{n_{0h}}{n_{0e}} = 1.0$, $\\frac{T_h}{T_c} = 101.695$")
p <- p +
        scale_x_continuous("$k_\\perp \\rho_h$", sec.axis = custom_secondary_axis, expand  = c(0,0)) +
        scale_y_continuous("$\\displaystyle \\frac{\\omega}{\\omega_{ce}}$", sec.axis = custom_secondary_axis) +
        scale_linetype_manual(name="Distribution", values=c("Maxwellian"="solid", "Double-Kappa"="dashed", "Cairns"="dotted", "VC"="dotdash")) +
        expand_limits(x = 0) +    # Make graph start at x = 0 (snug up with left margin)
        custom_theme

write_plot(p, filename, texOutput = T)
