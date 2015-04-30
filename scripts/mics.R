source("libgc.R")
OD.min = 0.01
OD.max = 0.5

gcs1 = read.gc("mics_plate1.txt")
names(gcs1)[-1] = paste0(names(gcs1)[-1],"_plate1")
gcs2 = read.gc("mics_plate2.txt")
names(gcs2)[-1] = paste0(names(gcs2)[-1],"_plate2")

#Join and normalize the data
cat("Joining and Normalizing...\n")
gcs = join(gcs1, gcs2)
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

area_control = mean(models$area[models$strain=="H2O"])
growth.plot = ggplot(models, aes(x=as.numeric(gsub("uM","",treatment)), 
    y=growth.rate, col=treatment, width=0.8)) +
    stat_smooth(method="lm", fill="slateblue", size=1, aes(group=1)) +
    geom_point(size=3, position=position_dodge(width=0.9)) +
	  scale_color_grey(start=0, end=0.6, name="") +
		xlab("concentration (μM)") + ylab("growth rate (1/h)") + pub_theme + 
    facet_wrap(~strain,nrow=2,scale="free_x") + theme(legend.position="none")

area.plot = ggplot(models, aes(x=as.numeric(gsub("uM","",treatment)), 
    y=area, color=treatment, width=0.8)) +
		geom_hline(yintercept=c(1,0.5)*area_control, linetype="dashed") +
    stat_smooth(method="lm", fill="slateblue", size=1, aes(group=1)) + 
    geom_point(size=3) + 
	  scale_color_grey(start=0, end=0.6, name="") +
		xlab("concentration (μM)") + ylab("area under curve (a.u.)") + pub_theme + facet_wrap(~strain,nrow=2, scale="free_x") +
		theme(legend.position="none")

# to get the original growth curves in nice formatting
names(gcs) = make.names(names(gcs), unique=T)
new = melt(gcs, id.vars="Time", variable.name="Probe", value.name="OD600")
fs = names_conv(new$Probe)
new  <- cbind(new, fs)
new$id = rep(1:(ncol(gcs)-1), each=nrow(gcs))

curves.plot = ggplot(new, aes(x=Time, y=OD600, col=treatment, group=id, shape=info)) + 
    geom_line(size=1, aes(linetype=info)) +
		#geom_line(stat="summary", fun.y=mean, size=1) + 
    facet_wrap(~strain) + 
		xlab("time (h)") + theme_bw()

cat(" - curves\n")		
ggsave(curves.plot, file="curves.svg", width=9, height=6)

cat(" - growth rates\n")
ggsave(growth.plot, file="growth.svg", width=9, height=4)

cat(" - area under curves\n")
ggsave(area.plot, file="area.svg", width=9, height=4)

capture.output(print(models), file="results.txt")
