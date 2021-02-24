installation.notes <- function()
{
	# shiver depends on python 2
	# brew install python@2 is no longer available, I did:	
	# pyenv install -v 2.7.17
	# pyenv global 3.8.2 2.7.17
	
	# pip3 install pyfastaq
	
	# pip install biopython	
	
	# brew install brewsci/bio/smalt
	
	# /Users/or105/.pyenv/versions/2.7.17/bin/python2.7 -m pip install biopython==1.76
	
	# mummer4 from
	# https://github.com/mummer4/mummer/blob/master/INSTALL.md
	
	# Trimmomatic
	# http://www.usadellab.org/cms/?page=trimmomatic
	# ~/sandbox/Trimmomatic-0.39/trimmomatic-0.39.jar
	
	# KMC 
	# https://github.com/refresh-bio/KMC/releases/tag/v3.1.1
}

one.sample.HIV <- function()
{	
	
	dir.shiver <- '/Users/or105/git/shiver'
	dir.project.base <- '/Users/or105/Box/OR_Work/2021/2021_HCV/tutorial_210115_'
	dir.init <- paste0(dir.project.base,'init_dir')
	dir.inputs <- paste0(dir.project.base,'settings')
	dir.contigs <- paste0(dir.project.base,'contigs')
	dir.fastq <- paste0(dir.project.base,'fastq')
	dir.out <- paste0(dir.project.base,'results')
	   
	
	#
	#	initialise project
	#
	
	#	initialise project - arguments
	file.config <- file.path('${DIR_INPUTS}','config.sh')
	file.primers <- file.path('${DIR_INPUTS}','primers_GallEtAl2012.fasta')
	file.adapters <- file.path('${DIR_INPUTS}','adapters_Illumina.fasta')
	file.refs <- file.path('${DIR_INPUTS}','HIV1_REF_2010_genome_DNA.fasta')
	
	#	initialise project - do it
	cmd <- ''
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")	
	cmd <- paste0(cmd,"DIR_INIT=",dir.init,"\n")	
	tmp <- paste0(
		file.path(dir.shiver, 'shiver_init.sh'),' ',
		'${DIR_INIT}',' ',
		file.config,' ',
		file.refs,' ',
		file.adapters,' ',
		file.primers)
	cmd <- paste0(cmd,tmp,"\n")
	cat(cmd)
	system(cmd)

	#
	#	assume IVA generated contigs
	#
	
	
	#
	#	process contigs
	#
	
	sample.id <- 'MysteryHIV'
	file.contigs <- file.path('${DIR_CONTIGS}','${SAMPLE_ID}_contigs.fasta')	
	dir.out.sample <- file.path(dir.out, '${SAMPLE_ID}')
	
	cmd <- ''
	cmd <- paste0(cmd,"CWD=$(pwd)\n")
	cmd <- paste0(cmd,"echo $CWD\n")
	cmd <- paste0(cmd,"SAMPLE_ID=",sample.id,"\n")
	cmd <- paste0(cmd,"DIR_INIT=",dir.init,"\n")
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")
	cmd <- paste0(cmd,"DIR_CONTIGS=",dir.contigs,"\n")
	cmd <- paste0(cmd,"mkdir -p ", dir.out.sample, "\n")
	cmd <- paste0(cmd,"cd ", dir.out.sample, "\n")
	tmp <- paste0(
			file.path(dir.shiver, 'shiver_align_contigs.sh'),' ',
			'${DIR_INIT}',' ',
			file.config,' ',
			file.contigs,' ',
			'${SAMPLE_ID}')
	cmd <- paste0(cmd,tmp,"\n")
	cmd <- paste0(cmd,"cd $CWD\n")
	cat(cmd)
	system(cmd)
	
	#
	#	continue mapping
	#
	file.blast <- '${SAMPLE_ID}.blast'
	file.processed.contigs <- '${SAMPLE_ID}_cut_wRefs.fasta'
	file.forward.fastq <- file.path('${DIR_FASTQ}', '${SAMPLE_ID}_1.fastq')
	file.backward.fastq <- file.path('${DIR_FASTQ}', '${SAMPLE_ID}_2.fastq')
			
	cmd <- ''
	cmd <- paste0(cmd,"CWD=$(pwd)\n")
	cmd <- paste0(cmd,"echo $CWD\n")	
	cmd <- paste0(cmd,"SAMPLE_ID=",sample.id,"\n")
	cmd <- paste0(cmd,"DIR_INIT=",dir.init,"\n")
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")
	cmd <- paste0(cmd,"DIR_FASTQ=",dir.fastq,"\n")
	cmd <- paste0(cmd,"DIR_CONTIGS=",dir.contigs,"\n")	
	cmd <- paste0(cmd,"cd ", dir.out.sample, "\n")
	tmp <- paste0(
			file.path(dir.shiver, 'shiver_map_reads.sh'),' ',
			'${DIR_INIT}',' ',
			file.config,' ',
			file.contigs,' ',
			'${SAMPLE_ID}', ' ',
			file.blast, ' ',
			file.processed.contigs, ' ',
			file.forward.fastq, ' ',
			file.backward.fastq
			)
	cmd <- paste0(cmd,tmp,"\n")	
	cmd <- paste0(cmd,"cd $CWD\n")
	cat(cmd)
	
	#	need trimmomatic
	#	script files will be easier to read if I define env variables with dir names
	'~/shiver/shiver_map_reads.sh MyInitDir config.sh contigs.fasta SID \
  SID.blast  SID_cut_wRefs.fasta  reads_1.fastq.gz  reads_2.fastq.gz
'
}

