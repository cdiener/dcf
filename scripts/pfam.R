DATA_DIR = "pfam"

get_prob = function(file_path)
{
	if( !file.exists( "pfam_data.txt" ) ) {
		output = system( sprintf("./predict ../examples 'CPP efficiency' %s/*.fasta > pfam_data.txt", 
						file_path), intern=TRUE )
	}
	
	df = read.table("pfam_data.txt", header=F)
	names(df) = c("family", "P_CPP", "P_eff", "P_all")

	return( df )
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
cat("Getting predictions...\n")
pfam = get_prob(DATA_DIR)
