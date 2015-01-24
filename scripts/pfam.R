DATA_DIR = "pfam/"

get_prob = function(file_path)
{
	output = system( sprintf("./predict %s ../examples CPP efficiency", 
					file_path), intern=TRUE )
	
	vals = strsplit(output, "\t")[[1]][1:3]
	
	return( vals )
}

get_pfam = function()
{
	setwd(DATA_DIR)
	
	# Download Pfam data if necessary
	if( length(dir(".", ".fasta$"))==0 && !file.exists("Pfam-A.fasta.gz") )
	{
		system("wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.fasta.gz")
	}
	
	# Unpack Pfam data if necessary
	if( file.exists("Pfam-A.fasta.gz") && !file.exists("Pfam-A.fasta") )
	{
		system("gzip -d Pfam-A.fasta.gz")
	}
	
	if( file.exists("Pfam-A.fasta") ) 
		system("../pfam_splitter Pfam-A.fasta && rm Pfam-A.fasta")
	
	setwd("..")
}

# The code that will be run

system("mkdir pfam 2> /dev/null")
get_pfam()

fasta_files = dir(DATA_DIR, ".fasta$")
families = sapply(fasta_files, sub, pattern="\\..+", replacement="")

probs = NULL

cat("Getting predictions\n")
for( f in fasta_files) {
	probs = rbind( probs, get_prob( paste0(DATA_DIR,f) ) )
	cat(".", file="")
}

pfam_data = data.frame(family=families, P.CPP=probs[,1], 
						P.eff=probs[,2], P.all=probs[,3])
						
save(pfam_data, "pfam_probs.Rd")

