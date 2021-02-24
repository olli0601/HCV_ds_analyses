make.adapters.primers.fasta <- function()
{
	file <- '~/Box/OR_Work/2021/2021_HCV/HCV1_210224_settings/HCV primers and index.csv'
	outdir <- dirname(file)
	
	dd <- read.csv(file, stringsAsFactors=FALSE)
	dd <- dd[dd$Seq!='',]	
	dd$Seq_Label <- paste0(gsub(' ','-',dd$Gene.Index),'_',gsub(' ','-',dd$Primer.Index))
	
	# select adapters and save as fasta
	adapters <- dd[grepl('adapter',dd$Gene.Index),]
	tmp <- rbind( paste0('>',adapters$Seq_Label, '\n'), paste0(adapters$Seq,'\n') ) 
	tmp <- paste0(as.vector(tmp),collapse='')
	cat(tmp, file=file.path(outdir, 'HCV_adapters.fasta'))
	
	# select primers and save as fasta
	primers <- dd[!grepl('adapter',dd$Gene.Index),]
	tmp <- rbind( paste0('>',primers$Seq_Label, '\n'), paste0(primers$Seq,'\n') ) 
	tmp <- paste0(as.vector(tmp),collapse='')
	cat(tmp, file=file.path(outdir, 'HCV_primers.fasta'))
}