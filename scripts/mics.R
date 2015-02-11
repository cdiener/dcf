source("libgc.R")
OD.min = 0.01
OD.max = 0.5

gcs = read.gc("MICs_24012015.txt")

#Normalize the data
cat("Normalizing...\n")
gcs = normalize(gcs)

# Get all model statistics
cat("Calculate Linear Regressions...\n")
models = lm.stats.multi(gcs)

# Diagnostic Plots
pdf("linear_parts.pdf", width=8, height=6)
plot(NULL, xlim=c(0,24), ylim=c(0.001,1.8), log="y", xlab="Time [h]", ylab="OD600")
 
gc = gcs[,-1]
dummy = apply(gc, 2, function(col) lines(gcs$Time, col, lwd=1, col="blue"))

abline(h=OD.min, lwd=2)
abline(h=OD.max, lwd=2)
dev.off()


# Rest of plots
cat("Plotting...\n")
source("plot_settings.R")
require(reshape2)

growth.plot = ggplot(models, aes(x=strain, y=growth.rate, fill=treatment, width=0.8)) +
		geom_bar(stat="summary", fun.y=mean, position=position_dodge(width=0.9)) +
		geom_errorbar(stat="summary", fun.data=mean_sdl, mult=1, width=0.5, position=position_dodge(width=0.9)) +  
		scale_fill_grey(name="") + scale_y_continuous(limits=c(0,0.45)) +  
		xlab("") + ylab("growth rate (1/h)") + pub_theme + 
		theme(axis.text.x = element_text(angle=45, hjust=1), legend.position=c(0.8,0.92), legend.direction="horizontal")

area.plot = ggplot(models, aes(x=strain, y=area, fill=treatment, width=0.8)) +
		geom_bar(stat="summary", fun.y=mean, position=position_dodge(width=0.9)) +
		geom_errorbar(stat="summary", fun.data=mean_sdl, mult=1, width=0.5, position=position_dodge(width=0.9)) +  
		scale_fill_grey(name="") + scale_y_continuous(limits=c(0,12)) +  
		xlab("") + ylab("area under curve (a.u.)") + pub_theme + 
		theme(axis.text.x = element_text(angle=45, hjust=1), legend.position=c(0.8,0.92), legend.direction="horizontal")

# to get the original growth curves in nice formatting
new = melt(gcs, id.vars="Time", variable.name="Probe", value.name="OD600")
fs = names_conv(new$Probe)
new  <- cbind(new, fs)

curves.plot = ggplot(new, aes(x=Time, y=OD600, col=treatment)) + 
		geom_linerange(stat="summary", fun.data=mean_sdl, mult=1, size=1, alpha=0.5 ) +
		geom_line(stat="summary", fun.y=mean, size=1) + facet_wrap(~strain) + 
		xlab("time (h)") + theme_bw()

cat(" - curves\n")		
ggsave(curves.plot, file="curves.svg", width=8, height=6)

cat(" - growth rates\n")
ggsave(growth.plot, file="growth.svg", width=9, height=5)

cat(" - area under curves\n")
ggsave(area.plot, file="area.svg", width=9, height=4)

capture.output(print(models), file="results.txt")
