N = 32
N_IT = 2000
AA = "ILWFVMYAPTSCGNDQEHKR"
AA = strsplit(AA,'')[[1]]
L = c(2,4,8,12,16)

before = list()
after = list()

bench = NULL

# Generate random sequences
system("mkdir seqs")

if( file.exists("bench_data.Rd") ){ 
	load("bench_data.Rd") 
} else {
	cat("Generating random sequences\n")
	for( i in 1:N )
	{
		pep = paste( sample(AA,16,replace=TRUE), collapse='' )
		write( pep, file=sprintf("seqs/random%d.txt", i) )
		system( sprintf("./modes seqs/random%d.txt '%d %d' %d ../examples CPP > /dev/null", 
				i, 1, 1, 1) )
		
		logf = read.table("log.txt", header=TRUE)
		bench = rbind( bench, c(0, 1-logf$best_E[1]) )
	}
	
	bench = as.data.frame(bench)
	names(bench)=c("length", "P_CPP")

	# Optimize random sequences
	for(l in 1:length(L))
	{
		cat( sprintf("\nTesting for linker size %d", L[l]) )

		for( i in 1:N )
		{
			system( sprintf("./modes seqs/random%d.txt '%d %d' %d ../examples CPP > /dev/null", 
					i, L[l], L[l], N_IT) )
			logf = read.table("log.txt", header=TRUE)
			bench = rbind( bench, c(L[l],1-logf$best_E[nrow(logf)]) )
			cat('.') 
		}
	}

	save(bench, file="bench_data.Rd")
}

source("plot_settings.R")
fragment_plot = ggplot(bench, aes(x=factor(length), y=P_CPP, group=length)) + 
				geom_boxplot() + xlab("fragment length (#aa)") + ylab("P(CPP)") +
				ylim(c(0,1)) + pub_theme

ggsave(fragment_plot, file="benchmark.svg", width=6, height=4)


system("./modes ../examples/alpha.txt '4 4' 500 ../examples CPP")
logf = read.table("log.txt", header=TRUE)

iter_plot = ggplot(logf, aes(x=iter, y=current_E)) + geom_line() + 
			xlab("iterations") + ylab("energy = 1-P(CPP)") + pub_theme
ggsave(iter_plot, file="sample_opt.svg", width=6, height=4)

pvals = sapply(L, function(l) t.test(bench$P_CPP[bench$length==0], 
				bench$P_CPP[bench$length==l], alternative="l")$p.value)
				
cat("Student t-test p-vals (alternative=less):\n")
print(pvals)
capture.output(pvals, file="t-test.txt")
