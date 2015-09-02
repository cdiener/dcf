DATA_DIR = "pfam"
PFAM_NAMESPACE = c(p = "http://pfam.xfam.org/")

get_prob = function(file_path)
{
    broken = !file.exists("pfam_data.txt")
    nf = as.numeric(system(sprintf("ls -1 %s/*.fasta | wc -l", file_path), intern=TRUE))
    if(file.exists("pfam_data.txt")) {
        nl = length(readLines("pfam_data.txt"))
        nf = as.numeric(system(sprintf("ls -1 %s/*.fasta | wc -l", file_path), intern=TRUE))
        broken = !(nl==nf)
        if(broken) {
            warning("pfam_data.txt is broken. Will be removed and regenerated!")
            file.remove("pfam_data.txt")
        }
    }
	if( broken ) {
		system( sprintf("./predict ../examples 'CPP efficiency' '%s'/*.fasta | 
            pv -s %d -l > pfam_data.txt", file_path, nf))
	}
	
	df = read.table("pfam_data.txt", header=F)
	names(df) = c("family", "P_CPP", "P_eff", "P_all")

	return( df )
}

get_pfam = function(release="current")
{
	setwd(DATA_DIR)
	
    if(release == "current") ver = "current_release"
    else ver = paste0("releases",release)
    
	# Download Pfam data if necessary
	if( length(dir(".", ".fasta$")) == 0 && !file.exists("Pfam-A.fasta.gz") )
	{
        write("Downloading PFAM...", file="")
		system(sprintf("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/%s/
            Pfam-A.fasta.gz", ver))
	}
	
	# Unpack Pfam data if necessary
	if( file.exists("Pfam-A.fasta.gz") && !file.exists("Pfam-A.fasta") )
	{
        write("Unpacking PFAM...", file="")
		system("gzip -d Pfam-A.fasta.gz")
	}
	
	if( file.exists("Pfam-A.fasta") ) 
    {
        write("Splitting PFAM...", file="")
		system("../pfam_splitter Pfam-A.fasta && rm Pfam-A.fasta")
    }
	
	setwd("..")
}

get_description = function(family_name)
{
	require(XML)
	x = xmlParse(sprintf("http://pfam.xfam.org/family/%s?output=xml",
			family_name) )
	desc = xmlChildren(x)[["pfam"]][["entry"]][["description"]]
	com = xmlChildren(x)[["pfam"]][["entry"]][["comment"]]
	
	com = gsub("\n", "", xmlValue(com))
	desc = gsub("\n", "", xmlValue(desc))
	
	return( paste(desc,com) )
}

get_go = function(family_name, type="function")
{
	require(XML)
	x = xmlParse(sprintf("http://pfam.xfam.org/family/%s?output=xml",
			family_name) )
	
	query = sprintf("//p:go_terms/p:category[@name='%s']/*", type)
	go_nodes = xpathApply(xmlRoot(x), query, namespaces=PFAM_NAMESPACE)
	if( is.null(go_nodes) ) return("N/A")
	
	return( paste(sapply(go_nodes, xmlValue), collapse=", ") )		
}

# The code that will be run

system("mkdir pfam 2> /dev/null")
get_pfam()
cat("Getting predictions. This will take a while...\n")
pfam = get_prob(DATA_DIR)
pfam = pfam[order(pfam$P_CPP, decreasing=TRUE),]
cpp = pfam[pfam$P_CPP>0.5,]

cat("Getting annotations\n")
desc = sapply(cpp$family, get_description)
cat('.')
go_funcs = sapply(cpp$family, get_go)
cat('.')
go_procs = sapply(cpp$family, get_go, type="process")
cat('.\n')
cpp = cbind(cpp, unlist(go_funcs), unlist(go_procs), unlist(desc))
names(cpp)[5:7] = c("GO_function", "GO_process", "description")

cat("Appearance of GO functions for CPP-like sequences\n:")
print( sort( table(go_funcs), decreasing=T ) )

cat("Appearance of GO processes for CPP-like sequences\n:")
print( sort( table(go_procs), decreasing=T ) )

write.csv(cpp, file="pfam_cpp.csv", row.names=F)

source("plot_settings.R")
pfam_hist = ggplot(pfam, aes(x=P_CPP)) + stat_bin(binwidth=0.05) + 
			geom_vline(xintercept=0.5, linetype="dashed") +
			scale_x_log10(breaks=10^(-4:0)) + xlab("P(CPP)") + pub_theme

ggsave(pfam_hist, file="pfam_hist.svg", width=6, height=4)
