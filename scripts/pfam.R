DATA_DIR = "pfam"

get_prob = function(file_path)
{
	output = system( sprintf("./predict %s ../examples CPP efficiency", file_path) )
	
	vals = strsplit(output, "\n")[[1]]
	vals = vals[length(vals)-1]
	vals = strsplit(vals, "\t")[[1]][1:3]
	
	return( vals )
}

install = function()
{
	write("Creating symbolic links...", file="")
	system( "ln -s ../pfam_splitter pfam_splitter" )
	system( "ln -s ../predict predict" )
}

get_pfam = function()
{
	setwd(DATA_DIR)
	
	# Download Pfam data if necessary
	if( !file.exists("Pfam-A.fasta.gz") && !file.exists("Pfam-A.fasta") )
	{
		system("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz")
	}
	
	# Unpack Pfam data if necessary
	if( file.exists("Pfam-A.fasta.gz") && !file.exists("Pfam-A.fasta") )
	{
		system("gzip e Pfam-A.fasta.gz")
	}
	
	if( file.exists("Pfam-A.fasta") ) 
		system("./pfam_splitter Pfam-A.fasta && rm Pfam-A.fasta")
	
	setwd("..")
}

# The code that will be run

install()
get_pfam()

fasta_files = dir(DATA_DIR, ".fasta$")
families = sapply(fasta_files, sub, patten="\\..+", replacement="")

probs = NULL

cat("Getting predictions", file="")
for( f in fasta_files) {
	probs = rbind(probs, get_probs(f))
	cat(".", file="")
}

pfam_data = data.frame(family=families, P.CPP=probs[,1], 
						P.eff=probs[,2], P.all=probs[,3])
						
save(pfam_data, "pfam_probs.Rd")

