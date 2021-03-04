May2018_v210224 <- function()
{	
	
	dir.shiver <- '/Users/or105/git/shiver'
	prog.trimmomatic <- '/Users/or105/sandbox/Trimmomatic-0.39/trimmomatic-0.39.jar'
	prog.iva <- 'iva'
	prog.docker.iva <- 'sangerpathogens/iva iva'
	prog.docker.trimmomatic <- '/Trimmomatic-0.38/trimmomatic-0.38.jar'
	smalt.k <- "15"
	smalt.s <- "3"
	smalt.y <- "0.7"
	regex.sample.id <- '(.*)_R([1-2])_[0-9]+.fastq.gz'
	regex.forward.id <- '(.*)_R1_[0-9]+.fastq.gz'
	regex.backward.id <- '(.*)_R2_[0-9]+.fastq.gz'
	
	dir.project.base <- '/Users/or105/Box/OR_Work/2021/2021_HCV/HCV1_210224_'
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
	file.primers <- file.path('${DIR_INPUTS}','HCV_primers.fasta')
	file.adapters <- file.path('${DIR_INPUTS}','HCV_adapters.fasta')
	file.refs <- file.path('${DIR_INPUTS}','HCV_genotypes_and_subtypes_v8-5-19.fasta')
	
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
	
	#
	#	read data to be processed
	#
	
	df <- data.table(F=list.files(dir.fastq))
	df[, SAMPLE_ID := gsub(regex.sample.id,'\\1',F)]
	df[, DIRECTION := NA_character_]
	tmp <- df[,which(grepl(regex.forward.id,F))]
	set(df, tmp, 'DIRECTION', 'fastq_fwd')
	tmp <- df[,which(grepl(regex.backward.id,F))]
	set(df, tmp, 'DIRECTION', 'fastq_bwd')
	stopifnot( df[, !any(is.na(DIRECTION))] )
	tmp <- df[, list(n=length(DIRECTION)), by='SAMPLE_ID']
	tmp2 <- subset(tmp, n!=2)[, unique(SAMPLE_ID)]
	if(length(tmp2))
		warning('Did not find bwd and fwd fastq for sample IDs. Excluding ', tmp2)
	df <- merge(df, subset(tmp, n==2, SAMPLE_ID), by='SAMPLE_ID')	
	df <- dcast.data.table(df, SAMPLE_ID~DIRECTION, value.var='F')
	
	
	#
	#	make contigs with IVA
	#
	
	data.dir <- "/Users/or105/Box/OR_Work/2021/2021_HCV/HCV1_210224_results/test"
	out.dir <- "docker_iva_out"
	cmd <- ''
	cmd <- paste0(cmd,"docker run --rm -it -v ",data.dir,":/data sangerpathogens/iva ")
	cmd <- paste0(cmd,"iva -f /data/reads_1.fq.gz -r /data/reads_2.fq.gz /data/",out.dir)
	
	i <- 1
		 
	cmd <- ''
	cmd <- paste0(cmd,"CWD=$(pwd)\n")
	cmd <- paste0(cmd,"echo $CWD\n")
	cmd <- paste0(cmd,"SAMPLE_ID=",df$SAMPLE_ID[i],"\n")
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")
	cmd <- paste0(cmd,"DIR_FASTQ=",dir.fastq,"\n")		
	cmd <- paste0(cmd,"DIR_CONTIGS=",dir.contigs,"\n")			
	cmd <- paste0(cmd,"DIR_CONTIGS_SAMPLE=${DIR_CONTIGS}/${SAMPLE_ID}\n")	
	cmd <- paste0(cmd,"FASTQ_FWD=",df$fastq_fwd[i],"\n")
	cmd <- paste0(cmd,"FASTQ_BWD=",df$fastq_bwd[i],"\n")
	cmd <- paste0(cmd,"mkdir -p ${DIR_CONTIGS_SAMPLE}\n")	
	cmd <- paste0(cmd,"cp -n ",file.adapters," ${DIR_CONTIGS_SAMPLE}/adapters.fasta\n")
	cmd <- paste0(cmd,"cp -n ",file.primers," ${DIR_CONTIGS_SAMPLE}/primers.fasta\n")	
	cmd <- paste0(cmd,"cp -n ${DIR_FASTQ}/${FASTQ_FWD} ${DIR_CONTIGS_SAMPLE}/${FASTQ_FWD}\n")
	cmd <- paste0(cmd,"cp -n ${DIR_FASTQ}/${FASTQ_BWD} ${DIR_CONTIGS_SAMPLE}/${FASTQ_BWD}\n")
	tmp <- paste0("docker run --rm -it",
			" -v ${DIR_CONTIGS_SAMPLE}:/data"," ",
			prog.docker.iva, " ",
			"--trimmomatic ",prog.docker.trimmomatic," ",
			"--adapters /data/adapters.fasta ",
			"--pcr_primers /data/primers.fasta ",
			"--smalt_k ",smalt.k," ",
			"--smalt_s ",smalt.s," ",
			"--smalt_id ",smalt.y," ",
			"-f /data/${FASTQ_FWD} ",
			"-r /data/${FASTQ_BWD} ",
			"/data/iva")
	cmd <- paste0(cmd,tmp,"\n")
	cmd <- paste0(cmd,"mv ${DIR_CONTIGS_SAMPLE}/iva/contigs.fasta ${DIR_CONTIGS}/${SAMPLE_ID}_contigs.fasta\n") 
	cmd <- paste0(cmd,"rm -rf ${DIR_CONTIGS_SAMPLE}\n")
	cmd <- paste0(cmd,"cd $CWD\n")
	cat(cmd)
	
	tmp <- paste0(prog.iva," ",
		"--trimmomatic ",prog.trimmomatic," ",
		"--adapters ",file.adapters," ",
		"--pcr_primers ",file.primers," ",
		"--smalt_k ",smalt.k," ",
		"--smalt_s ",smalt.s," ",
		"--smalt_id ",smalt.y," ",
		"-f ${DIR_FASTQ}/${FASTQ_FWD} ",
		"-r ${DIR_FASTQ}/${FASTQ_BWD} ",
		"${DIR_OUT}")
	cmd <- paste0(cmd,tmp,"\n")
	cmd <- paste0(cmd,"cd $CWD\n")
	cat(cmd)
	
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