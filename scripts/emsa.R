source("plot_settings.R")

normalize = function(exp_data)
{
	m = as.matrix(exp_data[,2:ncol(exp_data)])
	bg = t( replicate(nrow(m), m[1,]) )
	ref = t( replicate(nrow(m), m[nrow(m),]) )
	
	exp_data[,2:ncol(exp_data)] = 1 - (m-bg)/(ref-bg)
	
	return(exp_data[2:nrow(exp_data),])
}


emsa_1_1 = normalize( read.table("EMSA_1_1.txt", header=T) )
emsa_1_4 = normalize( read.table("EMSA_1_4.txt", header=T) )

require(reshape2)

emsa_1_1 = melt(emsa_1_1, id.vars="Name", variable.name="Repetition", value.name="intensity")   

emsa_1_4 = melt(emsa_1_4, id.vars="Name", variable.name="Repetition", value.name="intensity")   

emsa_1_1 = cbind(emsa_1_1, rep.int("1:1", nrow(emsa_1_1)) )
emsa_1_4 = cbind(emsa_1_4, rep.int("1:4", nrow(emsa_1_4)) )
names(emsa_1_1)[ncol(emsa_1_1)] = "charge_ratio"
names(emsa_1_4)[ncol(emsa_1_4)] = "charge_ratio"

emsa = rbind(emsa_1_1, emsa_1_4)

require(scales)

emsa_plot = ggplot(emsa, aes(x=Name, y=intensity, fill=charge_ratio, width=0.8)) + 
			geom_bar(stat="summary", fun.y=mean, position=position_dodge(width=0.9)) +
			geom_errorbar(stat="summary", fun.data=mean_sdl, mult=1, width=0.5, position=position_dodge(width=0.9)) +  
			scale_fill_grey(name="charge ratio") + scale_y_continuous(limits=c(0,1), breaks=seq(0,1,by=0.2), label=percent) +  
			xlab("") + ylab("DNA retention (%)") + pub_theme + scale_x_discrete(expand=c(0,1)) + 
			theme(axis.text.x = element_text(angle=45, hjust=1), legend.position=c(0.9,0.5))

ggsave(emsa_plot, file="emsa.svg", width=6, height=5)
