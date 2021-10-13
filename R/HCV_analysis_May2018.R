cmd.docker.iva <- function(sample.id,fastq.fwd,fastq.bwd,file.adapters,file.primers,prog.docker.iva,prog.docker.trimmomatic,dir.inputs,dir.fastq,dir.contigs,smalt.k,smalt.s,smalt.y)
{
	cmd <- ''
	cmd <- paste0(cmd,"CWD=$(pwd)\n")
	cmd <- paste0(cmd,"SMALT_K=",smalt.k,"\n")
	cmd <- paste0(cmd,"SMALT_S=",smalt.s,"\n")
	cmd <- paste0(cmd,"SMALT_ID=",smalt.y,"\n")
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")
	cmd <- paste0(cmd,"DIR_FASTQ=",dir.fastq,"\n")		
	cmd <- paste0(cmd,"DIR_CONTIGS=",dir.contigs,"\n")
	cmd <- paste0(cmd,"SAMPLE_ID=",sample.id,"\n")
	cmd <- paste0(cmd,"DIR_CONTIGS_SAMPLE=${DIR_CONTIGS}/${SAMPLE_ID}\n")	
	cmd <- paste0(cmd,"FASTQ_FWD=",fastq.fwd,"\n")
	cmd <- paste0(cmd,"FASTQ_BWD=",fastq.bwd,"\n")
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
			"--smalt_k ${SMALT_K} ",
			"--smalt_s ${SMALT_S} ",
			"--smalt_id ${SMALT_ID} ",
			"-f /data/${FASTQ_FWD} ",
			"-r /data/${FASTQ_BWD} ",
			"/data/iva")
	cmd <- paste0(cmd,tmp,"\n")
	cmd <- paste0(cmd,"mv ${DIR_CONTIGS_SAMPLE}/iva/contigs.fasta ${DIR_CONTIGS}/${SAMPLE_ID}_contigs.fasta\n") 
	cmd <- paste0(cmd,"rm -rf ${DIR_CONTIGS_SAMPLE}\n")
	cmd
}

cmd.shiver.align.contigs <- function(sample.id, contigs, file.config, prog.shiver.align, dir.init, dir.inputs, dir.contigs, dir.out)
{		
	cmd <- ''
	cmd <- paste0(cmd,"CWD=$(pwd)\n")
	cmd <- paste0(cmd,"echo $CWD\n")
	cmd <- paste0(cmd,"SAMPLE_ID=",sample.id,"\n")
	cmd <- paste0(cmd,"DIR_INIT=",dir.init,"\n")
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")
	cmd <- paste0(cmd,"DIR_CONTIGS=",dir.contigs,"\n")
	cmd <- paste0(cmd,"FILE_CONTIGS=",file.path('${DIR_CONTIGS}',contigs),"\n")
	cmd <- paste0(cmd,"DIR_OUT=",dir.out,"\n")
	cmd <- paste0(cmd,"DIR_OUT_SAMPLE=${DIR_OUT}/${SAMPLE_ID}\n")
	cmd <- paste0(cmd,"mkdir -p ${DIR_OUT_SAMPLE}\n")
	cmd <- paste0(cmd,"cd ${DIR_OUT_SAMPLE}\n")
	tmp <- paste0(
			prog.shiver.align,' ',
			'${DIR_INIT}',' ',
			file.config,' ',			
			'${FILE_CONTIGS} ',
			'${SAMPLE_ID}')
	cmd <- paste0(cmd,tmp,"\n")
	cmd <- paste0(cmd,"cd $CWD\n")
	cmd
}

cmd.shiver.map.contigs <- function(sample.id,fastq.fwd,fastq.bwd,contigs,ali.contigs,file.config,prog.shiver.map,dir.init, dir.inputs, dir.fastq,dir.contig,dir.out)
{
	cmd <- ''
	cmd <- paste0(cmd,"CWD=$(pwd)\n")
	cmd <- paste0(cmd,"echo $CWD\n")		
	cmd <- paste0(cmd,"DIR_INIT=",dir.init,"\n")
	cmd <- paste0(cmd,"DIR_INPUTS=",dir.inputs,"\n")
	cmd <- paste0(cmd,"DIR_FASTQ=",dir.fastq,"\n")
	cmd <- paste0(cmd,"DIR_CONTIGS=",dir.contigs,"\n")
	cmd <- paste0(cmd,"DIR_OUT=",dir.out,"\n")		
	cmd <- paste0(cmd,"SAMPLE_ID=",sample.id,"\n")
	cmd <- paste0(cmd,"FASTQ_FWD=",fastq.fwd,"\n")
	cmd <- paste0(cmd,"FASTQ_BWD=",fastq.bwd,"\n")
	cmd <- paste0(cmd,"FILE_CONTIGS=",contigs,"\n")
	cmd <- paste0(cmd,"FILE_ALI_CONTIGS=",ali.contigs,"\n")
	cmd <- paste0(cmd,"DIR_OUT_SAMPLE=${DIR_OUT}/${SAMPLE_ID}\n")
	cmd <- paste0(cmd,"cd ${DIR_OUT_SAMPLE}\n")
	tmp <- paste0(
			prog.shiver.map,' ',
			'${DIR_INIT}',' ',
			file.config,' ',
			'${DIR_CONTIGS}/${FILE_CONTIGS}',' ',
			'${SAMPLE_ID} ',
			'${SAMPLE_ID}.blast ',
			'${FILE_ALI_CONTIGS} ',
			'${DIR_FASTQ}/${FASTQ_FWD} ',
			'${DIR_FASTQ}/${FASTQ_BWD}'
	)		
	cmd <- paste0(cmd,tmp,"\n")	
	cmd <- paste0(cmd,"cd $CWD\n")
	cmd
}

May2018_v210224 <- function()
{	
	
	dir.shiver <- '/Users/or105/git/shiver'
	prog.trimmomatic <- '/Users/or105/sandbox/Trimmomatic-0.39/trimmomatic-0.39.jar'
	prog.iva <- 'iva'
	prog.docker.iva <- 'sangerpathogens/iva iva'
	prog.docker.trimmomatic <- '/Trimmomatic-0.38/trimmomatic-0.38.jar'
	prog.shiver.align <- file.path(dir.shiver, 'shiver_align_contigs.sh')
	prog.shiver.map <- file.path(dir.shiver, 'shiver_map_reads.sh')
	smalt.k <- "15"
	smalt.s <- "3"
	smalt.y <- "0.7"
	regex.sample.id <- '(.*)_R([1-2])_[0-9]+.fastq.gz'
	regex.forward.id <- '(.*)_R1_[0-9]+.fastq.gz'
	regex.backward.id <- '(.*)_R2_[0-9]+.fastq.gz'
	
	dir.project.base <- '/Users/or105/sandbox/HCV/HCV1_210224_'
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
	#	make contigs with IVA
	#	read data to be processed	
	#
	df <- data.table(F=list.files(dir.fastq, pattern='*.fastq.gz$'))
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
	#	build cmd	
	cmds <- vector('list', nrow(df))
	for(i in 1:nrow(df))
	{		
		cmd <- cmd.docker.iva(df$SAMPLE_ID[i],df$fastq_fwd[i],df$fastq_bwd[i],
				file.adapters,file.primers,
				prog.docker.iva,prog.docker.trimmomatic,
				dir.inputs,dir.fastq,dir.contigs,
				smalt.k,smalt.s,smalt.y)
		cmds[[i]] <- cmd
	}
	cmd <- paste0(cmds, collapse='\n#\n#      MAKE CONTIGS FOR SAMPLE\n#\n')
	cat(cmd, file=file.path(dir.out,'make_contigs.sh'))	
	system(paste0('chmod a+x ',file.path(dir.out,'make_contigs.sh')))
	
	#
	#	align IVA contigs against refs
	#	read data to be processed
	#	
	df <- data.table(F=list.files(dir.fastq, pattern='*.fastq.gz$'))
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
	#	read contigs
	tmp <- data.table(contigs=list.files(dir.contigs, pattern='*_contigs.fasta$'))
	tmp[, SAMPLE_ID := gsub('(.+)_contigs.fasta','\\1',contigs)]
	df <- merge(df, tmp, by='SAMPLE_ID')	
	#	build cmd
	cmds <- vector('list', nrow(df))
	for(i in 1:nrow(df))
	{
		cmd <- cmd.shiver.align.contigs(
			df$SAMPLE_ID[i], 
			df$contigs[i], 
			file.config, 
			prog.shiver.align, 
			dir.init, 
			dir.inputs, 
			dir.contigs, 
			dir.out)
		cmds[[i]] <- cmd
	}
	cmd <- paste0(cmds, collapse='\n#\n#      ALIGN CONTIGS FOR SAMPLE\n#\n')
	cat(cmd, file=file.path(dir.out,'align_contigs.sh'))
	system(paste0('chmod a+x ',file.path(dir.out,'align_contigs.sh')))
	
	#
	#	continue mapping
	#
	df <- data.table(F=list.files(dir.fastq, pattern='*.fastq.gz$'))
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
	#	read contigs
	tmp <- data.table(contigs=list.files(dir.contigs, pattern='*_contigs.fasta$'))
	tmp[, SAMPLE_ID := gsub('(.+)_contigs.fasta','\\1',contigs)]
	df <- merge(df, tmp, by='SAMPLE_ID')	
	#	read aligned contigs
	tmp <- data.table(ali_contigs_raw=list.files(dir.out, pattern='*_raw_wRefs.fasta$', recursive=TRUE))
	tmp[, SAMPLE_ID := dirname(ali_contigs_raw)]
	df <- merge(df, tmp, by='SAMPLE_ID', all.x=TRUE)	
	tmp <- data.table(ali_contigs_cut=list.files(dir.out, pattern='*_cut_wRefs.fasta$', recursive=TRUE))
	tmp[, SAMPLE_ID := dirname(ali_contigs_cut)]
	df <- merge(df, tmp, by='SAMPLE_ID', all.x=TRUE)	
	df[, ali_contigs_sel := ali_contigs_raw]
	tmp <- df[, which(!is.na(ali_contigs_cut))]
	set(df, tmp, 'ali_contigs_sel', df[tmp, ali_contigs_cut])
	df <- subset(df, !is.na(ali_contigs_sel))
	set(df, NULL, 'ali_contigs_sel', df[,basename(ali_contigs_sel)])
	#	build cmd
	cmds <- vector('list', nrow(df))
	for(i in 1:nrow(df))
	{
		cmd <- cmd.shiver.map.contigs(
			df$SAMPLE_ID[i],
			df$fastq_fwd[i],
			df$fastq_bwd[i],
			df$contigs[i],
			df$ali_contigs_sel[i],
			file.config,
			prog.shiver.map,
			dir.init, 
			dir.inputs, 
			dir.fastq,
			dir.contig,
			dir.out)
		cmds[[i]] <- cmd
	}
	cmd <- paste0(cmds, collapse='\n#\n#      MAP READS\n#\n')
	cat(cmd, file=file.path(dir.out,'map_reads.sh'))
	system(paste0('chmod a+x ',file.path(dir.out,'map_reads.sh')))	
}

inspect.fastq <- function()
{
	library(ShortRead)
	file <- '~/sandbox/HCV/HCV1_210224_results/010519-M10_S10_L001/010519-M10_S10_L001_R1_001.fastq'
	reads <- readFastq(file)
	sr <- sread(reads)
	sri <- id(reads)
	
	file2 <- '~/sandbox/HCV/HCV1_210224_results/010519-M10_S10_L001/010519-M10_S10_L001_R2_001.fastq'
	reads2 <- readFastq(file)
	sr2 <- sread(reads2)
	sri2 <- id(reads2)
}