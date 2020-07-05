#!/bin/bash
#########################################################
#  Reference-guided de novo assembly - ABYSS
# ====================================================
# Original by Heidi Lischer, 2015/2016
# Modified by Langqing Liu, 2020
#########################################################

# set variables #########################################
workPathFiles=/lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/SBSB
ref=/lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/refseq/SCEB/SCEB.fa
#refRed=/lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/refseq/SCEB/cut_SCEB_10K.fa
primerFile=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/ref-guided-assembly-pipeline/AdapterSeq_new.fa
primerFileMP=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/ref-guided-assembly-pipeline/AdapterSeqMP_new.fa

NThreads=30      # set the number of threads of every parallelizable step
maxReadLength=300
kmer=61         #define best K (need to be adapted)

# paired-end libraries -------------------
name=SBSB           # set name of your species
lib=(400)      # set insertion libraries
insLow=(150)     # lower bound of insertion size
insHigh=(550)  # upper bound of insertion size
libsd=(50)       # sd of insertion size

#readpath=/lustre/nobackup/WUR/ABGC/shared/Pig/Vasayan_Warty_pig_10X/Resequencing/runfile/tmpSVSV01M01
# list of files with forward reads according to lib array
reads1=(reads1.fq.gz)
# list of files with rewerse reads according to lib array
reads2=(reads2.fq.gz)
# short names of libraries
shortNames=(lib1) 

# set work path ---------------------------
workPath=${workPathFiles}/${name}_abyss

# log file
log=${workPath}/log_${name}_abyss.txt

# Programs --------------------------------
progFastQC=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/FastQC/fastqc
progTrimmomatic=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/Trimmomatic-0.32/trimmomatic-0.32.jar
progSamtools=/cm/shared/apps/samtools/gcc/64/1.5/bin/samtools
progVcfutils=/home/WUR/liu194/samtools/bcftools-1.9/bin/vcfutils.pl
progBcftools=/home/WUR/liu194/samtools/bcftools-1.9/bin/bcftools
progBamtools=/cm/shared/apps/bamtools/gcc/64/2.5.1/bin/bamtools
progBedtools=/cm/shared/apps/bedtools/gcc/64/2.26.0/bin/bedtools
progPicard=/cm/shared/apps/picard/2.7.1/picard.jar
progBowtie2=/cm/shared/apps/WUR/ABGC/bowtie2/bowtie2-2.2.3/bowtie2
progSeqtk=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/seqtk/seqtk
progAbyss=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/abyss-2.1.2/exe/abyss
progAmos=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/amos-3.1.0/bin
progNucmer=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/mummer-4.0.0beta2/nucmer
progGatk=/cm/shared/apps/WUR/ABGC/GATK/GATK3.5/GenomeAnalysisTK.jar
progSoapdenovo2=/cm/shared/apps/SHARED/SOAPdenovo2/SOAPdenovo2-bin-LINUX-generic-r240
progRemovShortSeq=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/RemoveShortSeq.jar
progGetBlocks=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/GetBlocks.jar
progFastaToAmos=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/FastaToAmos.jar
progWriteSoapConfig=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/WriteSoapConfig.jar
progFastaStats=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/FastaStats.jar
progSplitSeqLowCov=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/SplitSeqLowCov.jar
pregPrepREF=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/prepREF.sh
pregfasta_cutter=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/fasta_cutter
pregfastaFilterByGapLen=/lustre/nobackup/WUR/ABGC/liu194/analysis/Ref-Guideed-assemby/bin/fastaFilterByGapLen.pl
#########################################################



# run pipeline ##########################################

  #prepare reference
  cd ${workPathFiles}
  sh prepREF.sh ${pregfasta_cutter} ${pregfastaFilterByGapLen} ${ef}
  refRed=${workPathFiles}/cut_ref/fasta_cutter/results/cut_filtered_${ref}

# 1. Step: quality/adapter trimming and quality check:
#######################################################
  # quality check ----------
  echo "quality check of raw reads..."
  echo "quality check of raw reads..." > $log
  
  cd $workPathFiles
  fastqcOut=${workPathFiles}/FastQC_het
  mkdir -p ${fastqcOut}
    
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    ${progFastQC} -t ${NThreads} -o ${fastqcOut} ${reads1[i]} ${reads2[i]}
    echo "pass"
  done  
  
  
  # quality/adapter trimming --------
  # - remove Illumina adapters provided in the primer files
  # - remove leading and trailing low quality basses (<3) or N
  # - 4 base sliding window -> remove when average quality is < 15
  # - remove reads which are shorter than 40 bp
  echo "quality/adapter trimming..."
  echo "quality/adapter trimming..." >> $log
  
  trimOut=${workPathFiles}/Trim_het
  mkdir $trimOut
  
  read1TrimPair=()
  read1TrimUnPair=()
  read2TrimPair=()
  read2TrimUnPair=()
  for i in ${!lib[*]}  #for all indexes in the array
  do
    read1TrimPair[i]=${trimOut}/${shortNames[i]}_R1_trimPair.fastq 
    read1TrimUnPair[i]=${trimOut}/${shortNames[i]}_R1_trimUnPair.fastq 
    read2TrimPair[i]=${trimOut}/${shortNames[i]}_R2_trimPair.fastq 
    read2TrimUnPair[i]=${trimOut}/${shortNames[i]}_R2_trimUnPair.fastq 
    
    java -jar ${progTrimmomatic} PE -threads ${NThreads} ${reads1[i]} ${reads2[i]} ${read1TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimPair[i]} ${read2TrimUnPair[i]} ILLUMINACLIP:${primerFile}:2:30:7:5:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:40 >> $log
    echo "pass"
  done
  
    # quality check ----------
  echo "quality check of trimmed reads..."
  echo "quality check of trimmed reads..." >> $log
  
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    ${progFastQC} -t ${NThreads} -o ${fastqcOut} ${read1TrimPair[i]} ${read2TrimPair[i]} ${read1TrimUnPair[i]} ${read2TrimUnPair[i]}
    echo "pass"
  done
  
  
  

# 2. Step: map reads against reference 
#          and define blocks and superblocks
#######################################################
  # prepare Reference ---------
  echo "prepare reference..."
  echo "prepare reference..." >> $log
  
  #remove scaffolds shorter than 10 kb
  java -jar ${progRemovShortSeq} -i $ref -o ${refRed}.fa -length 10000 
  
  #create index files
  ${progSamtools} faidx ${refRed}.fa
  java -jar ${progPicard} CreateSequenceDictionary R=${refRed}.fa O=${refRed}.dict
  

  # map reads against reference ----------
  echo "run reference mapping..."
  echo "run reference mapping..." >> $log
  
  cd ${workPath}
  
  #merge unpaired files
  readTrimUnPair=${trimOut}/${name}_trimUnpair_mod.fastq
  cat ${read1TrimUnPair[*]} ${read2TrimUnPair[*]} > ${readTrimUnPair}
  cat ${trimOut}/*_R1_trimUnPair.fastq ${trimOut}/*_R1_trimUnPair.fastq > ${readTrimUnPair}
  libUnpair=Unpair
  
  # index reference file  
  echo "run reference mapping..."  
  ${progBowtie2}-build ${refRed}.fa ${refRed}
  
  mappedAll=()
  unmapped=()
  mapped=()
  mappedFiltered=()
  bowtieFailPair=()
  bowtieFailUnPair1=()
  bowtieFailUnPair2=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mappedAll[i]=${workPath}/${shortNames[i]}_all.sorted.bam
    unmapped[i]=${workPath}/${shortNames[i]}_unmapped.sorted.bam
    mapped[i]=${workPath}/${shortNames[i]}.sorted.bam
    mappedFiltered[i]=${workPath}/${shortNames[i]}.sorted.filtered.bam
    bowtieFailPair[i]=${workPath}/${shortNames[i]}_failPair.fastq
    bowtieFailUnPair1[i]=${workPath}/${shortNames[i]}_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${workPath}/${shortNames[i]}_failUnPairR2.fastq    
    (
      ${progBowtie2} --fast-local -p ${NThreads} -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${refRed} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -U ${read1TrimUnPair[i]},${read2TrimUnPair[i]} -S bowtie2.sam
      ${progSamtools} view -@ ${NThreads} -bS bowtie2.sam > bowtie2.bam
      ${progSamtools} sort -@ ${NThreads} -T ${shortNames[i]} -o ${mappedAll[i]} bowtie2.bam
      ${progSamtools} index -@ ${NThreads} ${mappedAll[i]}
    
      #filter unmapped reads    
      ${progSamtools} view -@ ${NThreads} -b -F 4 ${mappedAll[i]} > ${mapped[i]}
      ${progSamtools} index ${mappped[i]}
    
      #get unmapped reads
      ${progSamtools} view -@ ${NThreads} -b -f 4 ${mappedAll[i]} > ${unmapped[i]}
      ${progSamtools} view -@ ${NThreads} -b -f 9 ${unmapped[i]} > ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar ${progPicard} SamToFastq INPUT=${unmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${bowtieFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${bowtieFailPair[i]%.fastq}.2.fastq
      rm ${unmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -b -F 8 -f 64 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair1[i]}  
      ${progSamtools} view -b -F 8 -f 128 ${unmapped[i]} | ${progBamtools} convert -format fastq -out ${bowtieFailUnPair2[i]}
    
      ${progBamtools} stats -in ${mappedAll[i]} >> $log
      echo "--> ${mappedAll[i]}" >> $log
    
      #filter for mapping quality >=10    
      ${progSamtools} view -@ ${NThreads} -b -q 10 ${mapped[i]} > ${mappedFiltered[i]}
      ${progBamtools} stats -in ${mappedFiltered[i]} >> $log
      echo "--> ${mappedFiltered[i]}" >> $log
    
      #check insertion size
      java -jar ${progPicard} CollectInsertSizeMetrics R=${refRed}.fa I=${mapped[i]} O=${mapped[i]%.bam}_insertSize.txt HISTOGRAM_FILE=${mapped[i]%.bam}_insertSizeHist.pdf
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
    
  #merge alignment files
  mappedMerged=${workPath}/${name}_mate.sorted_mapped
  ${progSamtools} merge -@ ${NThreads} ${mappedMerged}.bam ${mapped[*]} ${mateMapped[*]}
  

  # get blocks and superblocks ----------
  echo "get blocks and superblocks..."
  echo "get blocks and superblocks..." >> $log
  
  #get coverage along genome
  covFile=${workPath}/${name}_mate_coverage.txt
  ${progBedtools} genomecov -ibam ${mappedMerged}.bam -bga > ${covFile}
  #only proparly paired reads
  ${progSamtools} view -bf 0x2 ${mappedMerged}.bam | ${progSamtools} sort - -n | ${progBedtools} bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | ${progBedtools} genomecov -i - -bga -g ${refRed}.fa.fai > ${covFile%.txt}Paired.txt
  
  #get blocks with a minimal coverage of 10 (paired-end) reads and create superblocks of at least 12000bp length and a minimal overlap of 300bp (max. overlap = 3*300bp)
  blocks=${workPath}/blocks.txt
  superblocks=${workPath}/superblocks.txt
  java -jar ${progGetBlocks} -i ${covFile} -paired ${covFile%.txt}Paired.txt -o ${blocks} -oSuper ${superblocks} -mCov 10 -sLength 12000 -sOverlap 300 -maxLength 100000

 for i in ${!lib[*]}
 do
 ${progSamtools} index -@ ${NThreads} ${mapped[i]}
 done 

# 3. Step: do deNovo assembly within superblocks
#######################################################
  echo "deNovo assembly within superblocks..."
  echo "deNovo assembly within superblocks..." >> $log
  blocks=${workPath}/blocks.txt
  superblocks=${workPath}/superblocks.txt
  cd ${workPath}
  mappedAll=()
  unmapped=()
  mapped=()
  mappedFiltered=()
  bowtieFailPair=()
  bowtieFailUnPair1=()
  bowtieFailUnPair2=()
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mappedAll[i]=${workPath}/${shortNames[i]}_all.sorted.bam
    unmapped[i]=${workPath}/${shortNames[i]}_unmapped.sorted.bam
    mapped[i]=${workPath}/${shortNames[i]}.sorted.bam
    mappedFiltered[i]=${workPath}/${shortNames[i]}.sorted.filtered.bam
    bowtieFailPair[i]=${workPath}/${shortNames[i]}_failPair.fastq
    bowtieFailUnPair1[i]=${workPath}/${shortNames[i]}_failUnPairR1.fastq
    bowtieFailUnPair2[i]=${workPath}/${shortNames[i]}_failUnPairR2.fastq 
  done
  abyssRes=${workPath}/abyssResults
  mkdir ${abyssRes}
  mkdir ${abyssRes}/contigs
  mkdir ${abyssRes}/scaffolds
  
  #split superbolcks.txt into ${NThreads} files to run it in parallel
  size=$(($(wc -l < ${superblocks})/$((NThreads))+1))
  split -l $size -d ${superblocks} ${superblocks%.txt}_

    for (( j=10; j<$((NThreads)); j++ ))
  do
    mv ${superblocks%.txt}_${j} ${superblocks%.txt}_0${j}
  done

  count=0
  for (( j=0; j<$((NThreads)); j++ ))
  do
    file=${superblocks%.txt}_0${j}
    fileout=${superblocks%.txt}_0${j}_run.sh
    logout=${superblocks%.txt}_0${j}.log
    
    array=(${file//_/ })
    number=${array[${#array[*]}-1]}
    if [ ${number} != "00" ]
    then
      number=`echo $number|sed 's/^0*//'`
    fi
    
    if [[ $number =~ ^[0-9]+$ ]]
    then
      blockNb=$(($number*$size+1))
    else
      blockNb=1
    fi
    
    printf "#"'!'"/bin/bash\n"                 > ${fileout}
    printf "\n"                                >> ${fileout}
    printf "mkdir ${file}_temp\n"              >> ${fileout}
    printf "cd ${file}_temp\n"                 >> ${fileout}
    printf "\n" >> ${fileout}

    printf "blockNb=$blockNb\n" >> ${fileout}
    printf "start=\`date +%%s\`\n" >> ${fileout}
    printf "for block in \$(cat ${file})\n" >> ${fileout}
    printf "do\n" >> ${fileout}
    printf "  echo \$blockNb >> $logout\n" >> ${fileout}
    printf "\n" >> ${fileout}
    printf "  #extract sequence names within specified region\n" >> ${fileout}
    
    seqNames=()
    seqBam=()
    subSeq1=()
    subSeq2=()
    for i in ${!lib[*]}
    do
      seqNames[i]=sequences_${shortNames[i]}
      seqBam[i]=sequences_${shortNames[i]}.bam
      subSeq1[i]=subseq_${shortNames[i]}_R1.fastq
      subSeq2[i]=subseq_${shortNames[i]}_R2.fastq
      printf "  ${progSamtools} view -b ${mapped[i]} \$block | ${progSamtools} sort - -T ${shortNames[i]} -no ${seqBam[i]}\n" >> ${fileout}
      printf "  ${progBedtools} bamtofastq -i ${seqBam[i]} -fq 1_${subSeq1[i]} -fq2 1_${subSeq2[i]}\n" >> ${fileout}

      #extract paired reads with one pair unmapped
      printf "  ${progSamtools} view  -b -f 72 ${seqBam[i]} | ${progBamtools} convert -format fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq1[i]}\n" >> ${fileout}
      printf "  ${progSamtools} view  -f 72 ${seqBam[i]} | cut -f 1 | awk \'{print \$0\"/2\"}\' > ${seqNames[i]}_R2.txt\n" >> ${fileout}
      printf "  ${progSeqtk} subseq ${bowtieFailUnPair2[i]} ${seqNames[i]}_R2.txt | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 2_${subSeq2[i]}\n" >> ${fileout}
  
      printf "  ${progSamtools} view  -b -f 136 ${seqBam[i]} | ${progBamtools} convert -format fastq | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq1[i]}\n" >> ${fileout}
      printf "  ${progSamtools} view  -f 136 ${seqBam[i]} | cut -f 1 | awk \'{print \$0\"/1\"}\' > ${seqNames[i]}_R1.txt\n" >> ${fileout}
      printf "  ${progSeqtk} subseq ${bowtieFailUnPair1[i]} ${seqNames[i]}_R1.txt | paste - - - - | sort -k1,1 -t \" \" | tr \"\\\t\" \"\\\n\" > 3_${subSeq2[i]}\n" >> ${fileout}
  
      printf "  cat 1_${subSeq1[i]} 2_${subSeq1[i]} 3_${subSeq1[i]} > ${subSeq1[i]}\n" >> ${fileout}
      printf "  cat 1_${subSeq2[i]} 2_${subSeq2[i]} 3_${subSeq2[i]} > ${subSeq2[i]}\n" >> ${fileout}
      printf "\n" >> ${fileout}
    done

    printf "  #deNovo assembly (k41,k51,k61,k71,k81)------\n" >> ${fileout}
    printf "  #ABySS\n" >> ${fileout}
    for i in ${!lib[*]}
    do
      if [ $i == 0 ]
      then
        abyssLib=pe${lib[i]}
        abyssPaired="pe${lib[i]}=\"${subSeq1[i]} ${subSeq2[i]}\""
      else
        abyssLib="${abyssLib} pe${lib[i]}"
        abyssPaired="${abyssPaired} pe${lib[i]}=\"${subSeq1[i]} ${subSeq2[i]}\""
      fi      
    done
    #printf "  unset NSLOTS LSB_DJOB_NUMPROC SLURM_NTASKS\n" >> ${fileout}
    printf "  ${progAbyss}-pe k=41 name=sblock\${blockNb}_41 lib=\"${abyssLib}\" ${abyssPaired}\n" >> ${fileout}
    printf "  ${progAbyss}-pe k=51 name=sblock\${blockNb}_51 lib=\"${abyssLib}\" ${abyssPaired}\n" >> ${fileout}
    printf "  ${progAbyss}-pe k=61 name=sblock\${blockNb}_61 lib=\"${abyssLib}\" ${abyssPaired}\n" >> ${fileout}
    printf "  ${progAbyss}-pe k=71 name=sblock\${blockNb}_71 lib=\"${abyssLib}\" ${abyssPaired}\n" >> ${fileout}
    printf "  ${progAbyss}-pe k=81 name=sblock\${blockNb}_81 lib=\"${abyssLib}\" ${abyssPaired}\n" >> ${fileout}
    printf "  if [ ! -f sblock\${blockNb}_61-contigs.fa ]\n" >> ${fileout}
    printf "  then\n" >> ${fileout}
    printf "    echo \"\$blockNb abyss failed\" >> $logout\n" >> ${fileout}
    printf "  fi\n" >> ${fileout}
    printf "  cp *-contigs.fa ${abyssRes}/contigs/.\n" >> ${fileout}
    printf "  cp *-scaffolds.fa ${abyssRes}/scaffolds/.\n" >> ${fileout}
    printf "\n" >> ${fileout}
           
    printf "  rm -rf ${file}_temp/*\n" >> ${fileout}
    printf "  ((blockNb++))\n" >> ${fileout}
    printf "done\n" >> ${fileout}
    printf "\n" >> ${fileout}   
    printf "rm -rf ${file}_temp\n" >> ${fileout}
    printf "end=\`date +%%s\`\n" >> ${fileout}
    printf "echo \$((end-start))\n" >> ${fileout}
    printf "\n" >> ${fileout}
    
    chmod +x ${fileout}
    ${fileout} &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait


  # deNovo assembly of unassembled reads ----------
  echo "deNovo assembly of unassembled reads..."
  echo "deNovo assembly of unassembled reads..." >> $log
    
  unassFolder=${workPath}/Unassembled
  mkdir ${unassFolder}
  cd ${unassFolder}
  
  #merge unassembled unpaired files
  bowtieFailUnpairMerged=${workPath}/${name}_failUnp.fastq
  cat ${workPath}/*_failUnPairR1.fastq ${workPath}/*_failUnPairR2.fastq > $bowtieFailUnpairMerged
  
  for i in ${!lib[*]}
  do
    if [ $i == 0 ]
    then
      abyssLib=pe${lib[i]}
      abyssPaired="pe${lib[i]}=\"${bowtieFailPair[i]%.fastq}.1.fastq ${bowtieFailPair[i]%.fastq}.2.fastq\""
    else
      abyssLib="${abyssLib} pe${lib[i]}"
      abyssPaired="${abyssPaired} pe${lib[i]}=\"${bowtieFailPair[i]%.fastq}.1.fastq ${bowtieFailPair[i]%.fastq}.2.fastq\""
    fi      
  done
  
  count=0
  for (( j=41; j<=81; j+=10 ))
  do
    mkdir ${unassFolder}/k${j}
    cd ${unassFolder}/k${j}
    printf "#"'!'"/bin/bash\n" > runAbyss_k${j}.sh
    # np= CPU / 5
    printf "${progAbyss}-pe np=4 k=${j} name=Unass_${j} lib=\"${abyssLib}\" ${abyssPaired} se=${bowtieFailUnpairMerged}" >> runAbyss_k${j}.sh
    chmod +x runAbyss_k${j}.sh
    ./runAbyss_k${j}.sh &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
  

# 4. Step: get non-redundant supercontigs
####################################################### 
  echo "get supercontigs..."
  echo "get supercontigs..." >> $log
  cd ${workPath}
  
  #merge deNovo assembled superblocks into one FASTA file
  #list is too long for directly cat
  abyssContigs=${abyssRes}/abyssContigs.fa
  for (( j=41; j<=81; j+=10 ))
  do
    rm ${abyssRes}/${j}_abyssContigs.fa
    find ${abyssRes}/contigs/ -name "*${j}-contigs.fa" | while read F; do cat ${F} >> ${abyssRes}/${j}_abyssContigs.fa; done
    #ls ${abyssRes}/contigs/*${j}-contigs.fa | xargs -n 30 -P 6 cat > ${abyssRes}/${j}_abyssContigs.fa
    #cat ${abyssRes}/contigs/*${j}-contigs.fa > ${abyssRes}/${j}_abyssContigs.fa
  done   
  cat ${abyssRes}/*_abyssContigs.fa > ${abyssContigs}
  
  #merge Unassembled files
  orgUnass=${unassFolder}/Unassembled_Abyss.fa
  cat ${unassFolder}/k*/*-contigs.fa > ${orgUnass}

  #remove short seq (<500)
  echo "remove seq < 500 in ${orgUnass}" >> $log
  unass500=${abyssUnass%.fa}_500.fa
  java -jar ${progRemovShortSeq} -i ${orgUnass} -o ${unass500} -length 500  >> $log
  
  #merge all files
  superblockSeq=${workPath}/deNovo_Superblocks.fa
  cat ${abyssContigs} ${unass500} > ${superblockSeq}
    
  #remove short seq (<200)
  echo "remove seq < 200 in ${superblockSeq}" >> $log
  superblockSeq200=${superblockSeq%.fa}_200.fa
  java -jar ${progRemovShortSeq} -i ${superblockSeq} -o ${superblockSeq200} -length 200 -u  >> $log

  #remove redundency with AMOScmp
  amosFolder=${workPath}/AMOScmp
  mkdir ${amosFolder}
  cd ${amosFolder}

  #assemble all assembled superblocks with AMOScmp to supercontigs (with the help of reference)
  #changed parameters in AMOScmp: (casm-layout -t 1000 (maximum ignorable trim length), make-consensus -o 10 (minimum overlap base)) 
  superblockSeqAmos=${superblockSeq200%.fa}_Amos.afg 
  java -jar ${progFastaToAmos} -i ${superblockSeq200} -o ${superblockSeqAmos}
  supercontigs=Amos_supercontigs
  amosSupercontigs=${amosFolder}/${supercontigs}.fasta
  
  echo "run AMPScmp..." >> $log
  #${progAmos}AMOScmp -D TGT=${superblockSeqAmos} -D REF=${refRed}.fa ${supercontigs}
  
  # running AMPScmp step by step and use multithread nucmer to spead it up
  ## Building AMOS bank
  echo "  build AMPS bank..." >> $log
  ${progAmos}/bank-transact -c -z -b ${supercontigs}.bnk -m ${superblockSeqAmos}

  ## Collecting clear range sequences
  echo "  clear range sequences..." >> $log
  ${progAmos}/dumpreads ${supercontigs}.bnk > ${supercontigs}.seq

  ## Running nucmer
  echo "  run nucmer..." >> $log
  #${progNucmer} --maxmatch --threads=${NThreads} --prefix=${supercontigs} ${refRed}.fa ${supercontigs}.seq
  rm Amos_supercontigs.delta
  ${progNucmer} --threads=${NThreads} --prefix=${supercontigs} ${refRed}.fa ${supercontigs}.seq

  ## Running layout
  echo "  run layout..." >> $log
  ${progAmos}/casm-layout -t 1000 -U ${supercontigs}.layout -C ${supercontigs}.conflict -b ${supercontigs}.bnk ${supercontigs}.delta

  ## Running consensus
  echo "  run consensus..." >> $log
  ${progAmos}/make-consensus -o 10 -B -b ${supercontigs}.bnk

  ## Outputting contigs
  echo "  output contigs..." >> $log
  ${progAmos}/bank2contig ${supercontigs}.bnk > ${supercontigs}.contig

  ## Outputting fasta
  echo "  output fasta..." >> $log
  ${progAmos}/bank2fasta -b ${supercontigs}.bnk > ${supercontigs}.fasta
      
 
 
# 5. Step: map reads on supercontigs
#          and de novo assemble unmapped reads
####################################################### 
  echo "map reads on supercontigs and correct them..."
  echo "map reads on supercontigs and correct them..." >> $log
  
  #make seqnames unique
  amosSupercontigsUnique=${amosSupercontigs%.fasta}_unique.fa
  java -jar ${progRemovShortSeq} -i ${amosSupercontigs} -o ${amosSupercontigsUnique} -length 1 -u
   
  #get statistics
  echo ${amosSupercontigsUnique} >> $log
  java -jar ${progFastaStats} -i ${amosSupercontigsUnique} -min 200 >> $log
  
  #prepare reference
  ${progBowtie2}-build ${amosSupercontigsUnique} ${amosSupercontigsUnique%.fa}
  
  supercontMappedAll=()
  supercontUnmapped=()
  supercontFailPair=()
  supercontFailUnpair=()
  supercontMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    supercontMappedAll[i]=${amosFolder}/${shortNames[i]}_all.sorted.bam
    supercontUnmapped[i]=${amosFolder}/${shortNames[i]}_unmapped.sorted.bam
    supercontFailPair[i]=${amosFolder}/${shortNames[i]}_failPair.fastq
    supercontFailUnpair[i]=${amosFolder}/${shortNames[i]}_failUnp.fastq
    supercontMappedFiltered[i]=${amosFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p ${NThreads} -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${amosSupercontigsUnique%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -S bowtie2.sam
      ${progSamtools} view -@ ${NThreads} -bS bowtie2.sam > bowtie2.bam
      ${progSamtools} sort -@ ${NThreads} -T ${shortNames[i]} -o ${supercontMappedAll[i]} bowtie2.bam
      ${progSamtools} index -@ ${NThreads} ${supercontMappedAll[i]}
      #get unmapped reads
      ${progSamtools} view -@ ${NThreads} -b -f 4 ${supercontMappedAll[i]} > ${supercontUnmapped[i]}
      ${progSamtools} view -b -@ ${NThreads} -f 9 ${supercontUnmapped[i]} > ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      java -jar ${progPicard} SamToFastq INPUT=${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam FASTQ=${supercontFailPair[i]%.fastq}.1.fastq SECOND_END_FASTQ=${supercontFailPair[i]%.fastq}.2.fastq
      rm ${supercontUnmapped[i]%.sorted.bam}_pair.sorted.bam
      ${progSamtools} view -b -F 8 ${supercontUnmapped[i]} | ${progBamtools} convert -format fastq -out ${supercontFailUnpair[i]}
    
      ${progBamtools} stats -in ${supercontMappedAll[i]} >> $log
      echo "--> ${supercontMappedAll[i]}" >> $log
     
      #filter for mapping quality >=10    
      ${progSamtools} view -@ ${NThreads} -b -F 4 -q 10 ${supercontMappedAll[i]} > ${supercontMappedFiltered[i]}
      ${progBamtools} stats -in ${supercontMappedFiltered[i]} >> $log
      echo "--> ${supercontMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
       

  # deNovo assemble unassembled reads ----------
  echo "deNovo assemble unassembled reads..."
  echo "deNovo assemble unassembled reads..." >> $log
  
  supercontFailUnpairMerged=${amosFolder}/${name}_failUnp.fastq
  cat ${amosFolder}/*_failUnp.fastq > ${supercontFailUnpairMerged}
  
  supercontUnassFolder=${amosFolder}/Unassembled
  mkdir ${supercontUnassFolder}
  cd ${supercontUnassFolder}
  
  #deNovo assembly (k41,k51,k61,k71,k81)------
  #ABySS
  for i in ${!lib[*]}
  do
    if [ $i == 0 ]
    then
      abyssLib=pe${lib[i]}
      abyssPaired="pe${lib[i]}=\"${supercontFailPair[i]%.fastq}.1.fastq ${supercontFailPair[i]%.fastq}.2.fastq\""
    else
      abyssLib="${abyssLib} pe${lib[i]}"
      abyssPaired="${abyssPaired} pe${lib[i]}=\"${supercontFailPair[i]%.fastq}.1.fastq ${supercontFailPair[i]%.fastq}.2.fastq\""
    fi      
  done
  
  mkdir ${supercontUnassFolder}/k${kmer}
  cd ${supercontUnassFolder}/k${kmer}
  printf "#"'!'"/bin/bash\n" > runAbyss_k${kmer}.sh
  #printf "source /lustre/nobackup/WUR/ABGC/liu194/analysis/bin/miniconda3/bin/activate /lustre/nobackup/WUR/ABGC/liu194/analysis/bin/miniconda3/envs/ABYSS\n" >> runAbyss_k${kmer}.sh
  printf "${progAbyss}-pe k=${kmer} name=Unass_${kmer} lib=\"${abyssLib}\" ${abyssPaired} se=${supercontFailUnpairMerged}\n" >> runAbyss_k${kmer}.sh
  #printf "source  /lustre/nobackup/WUR/ABGC/liu194/analysis/bin/miniconda3/bin/deactivate\n" >> runAbyss_k${kmer}.sh
  chmod +x runAbyss_k${kmer}.sh
  ./runAbyss_k${kmer}.sh
  cd ${supercontUnassFolder}
BLOCK
  #remove contigs shorter than 200 bp
  cd ${supercontUnassFolder}
  supercontSeqUnass=${supercontUnassFolder}/Unass-contigs_200.fa
  java -jar ${progRemovShortSeq} -i k${kmer}/Unass_${kmer}-contigs.fa -o ${supercontSeqUnass} -length 100 >> $log
  echo ${supercontSeqUnass} >> $log
  java -jar ${progFastaStats} -i ${supercontSeqUnass} -min 200 >> $log
    
 
 
# 6. Step: map reads to all supercontics and correct them 
########################################################## 
  echo "merge contigs..." 
  echo "merge contigs..." >> $log
  
  mergedFolder=${workPath}/merged_corr
  mkdir $mergedFolder
  cd $mergedFolder
  
  merged=${mergedFolder}/${name}_supercontSeq_Unass.fa
  cat ${amosSupercontigsUnique} ${supercontSeqUnass} > $merged
  
  echo "${merged}" >> $log
  java -jar ${progFastaStats} -i ${merged} -min 200 >> $log
  
    ${progBowtie2}-build ${merged} ${merged%.fa}
  
  mergedMappedAll=()
  mergedMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mergedMappedAll[i]=${mergedFolder}/${shortNames[i]}_all.sorted.bam
    mergedMappedFiltered[i]=${mergedFolder}/${shortNames[i]}.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p ${NThreads} -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${merged%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -S bowtie2.sam
      ${progSamtools} view -@ ${NThreads} -bS bowtie2.sam > bowtie2.bam
      ${progSamtools} sort -@ ${NThreads} -T ${shortNames[i]} -o ${mergedMappedAll[i]} bowtie2.bam
      ${progSamtools} index -@ ${NThreads} ${mergedMappedAll[i]}

      ${progBamtools} stats -in ${mergedMappedAll[i]} >> $log
      echo "--> ${mergedMappedAll[i]}" >> $log
     
      #filter for mapping quality >=10
      ${progSamtools} view -@ ${NThreads} -b -F 4 -q 10 ${mergedMappedAll[i]} > ${mergedMappedFiltered[i]}
    
      ${progBamtools} stats -in ${mergedMappedFiltered[i]} >> $log
      echo "--> ${mergedMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait
  
  
  # error correction ----------  
  # add RG header of second file (lost while merging)
  for i in ${!lib[*]}  #for all indexes in the array
  do
    echo -e "@RG\tID:${shortNames[i]}.filtered.sorted\tPL:illumina\tPU:${lib[i]}\tLB:${lib[i]}\tSM:${shortNames[i]}" >> rg
  done
  ${progSamtools} view -H ${mergedMappedFiltered[1]} | cat - rg > header
      
  #realign reads
  mergedMappedMerged=${mergedFolder}/${name}.filtered_RG.sorted.bam
  ${progSamtools} merge -@ ${NThreads} -r -h header ${mergedMappedMerged} ${mergedMappedFiltered[*]}
  rm rg header
 
  ${progSamtools} index -@ ${NThreads} ${mergedMappedMerged}
  ${progSamtools} faidx ${merged}
  java -jar ${progPicard} CreateSequenceDictionary R=${merged} O=${merged%.fa}.dict
  java -jar ${progGatk} -nt ${NThreads} -T RealignerTargetCreator -R ${merged} -I ${mergedMappedMerged} -o target_intervals.list
  mergedMappedMergedReal=${mergedMappedMerged%.bam}_realigned.bam
  java -jar ${progGatk} -T IndelRealigner -R ${merged} -I ${mergedMappedMerged} -targetIntervals target_intervals.list -o ${mergedMappedMergedReal}

  #get alternative seq
  mergedCorr=${merged%.fa}_corr.fq
  ${progSamtools} mpileup -uf ${merged} ${mergedMappedMergedReal} | ${progBcftools} call -c - | ${progVcfutils} vcf2fq -d 1 > ${mergedCorr}    
  
  #remove start and end N
  mergedCorrWN=${mergedCorr%.fq}WN.fa 
  echo ${mergedCorrWN} >> $log
  java -jar ${progRemovShortSeq} -i ${mergedCorr} -o ${mergedCorrWN} -length 100 -n -fq >> $log
  
  #get statistics
  echo ${mergedCorrWN} >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN} -min 200 >> $log
  

  # split sequences at places with no coverage ----------
  ${progBowtie2}-build ${mergedCorrWN} ${mergedCorrWN%.fa}
  
  mergedCorrMappedAll=()
  mergedCorrMappedFiltered=()
  count=0
  for i in ${!lib[*]}  #for all indexes in the array
  do 
    mergedCorrMappedAll[i]=${mergedFolder}/${shortNames[i]}_corrWN_all.sorted.bam
    mergedCorrMappedFiltered[i]=${mergedFolder}/${shortNames[i]}_corrWN.filtered.sorted.bam
    (
      ${progBowtie2} --sensitive -p ${NThreads} -q --phred33 -I ${insLow[i]} -X ${insHigh[i]} -x ${mergedCorrWN%.fa} -1 ${read1TrimPair[i]} -2 ${read2TrimPair[i]} -S bowtie2.sam
      ${progSamtools} view -@ ${NThreads} -bS bowtie2.sam > bowtie2.bam
      ${progSamtools} sort -@ ${NThreads} -T ${shortNames[i]} -o ${mergedCorrMappedAll[i]} bowtie2.bam
      ${progSamtools} index -@ ${NThreads} ${mergedCorrMappedAll[i]}

      ${progBamtools} stats -in ${mergedCorrMappedAll[i]} >> $log
      echo "--> ${mergedCorrMappedAll[i]}" >> $log
    
      #filter for mapping quality >=10
      ${progSamtools} view -@ ${NThreads} -b -F 4 -q 10 ${mergedCorrMappedAll[i]} > ${mergedCorrMappedFiltered[i]}
    
      ${progBamtools} stats -in ${mergedCorrMappedFiltered[i]} >> $log
      echo "--> ${mergedCorrMappedFiltered[i]}" >> $log
    ) &
    let count+=1
    [[ $((count%${NThreads})) -eq 0 ]] && wait
  done
  wait

  mergedCorrMappedFilteredMerged=${mergedFolder}/${name}_corrWN.filtered.sorted.bam
  #${progSamtools} merge ${mergedCorrMappedFilteredMerged} ${mergedCorrMappedFiltered[*]}
  cp ${mergedCorrMappedFiltered[*]} ${mergedCorrMappedFilteredMerged} 
  ${progBedtools} genomecov -ibam ${mergedCorrMappedFilteredMerged} -bga > ${mergedCorrWN%.fa}_filteredCov.txt
  #only proparly paired reads
  ${progSamtools} faidx $mergedCorrWN
  ${progSamtools} view -bf 0x2 ${mergedCorrMappedFilteredMerged} | ${progSamtools} sort - -n | ${progBedtools} bamtobed -i - -bedpe | awk '$1 == $4' | cut -f 1,2,6 | sort -k 1,1 | ${progBedtools} genomecov -i - -bga -g ${mergedCorrWN}.fai  > ${mergedCorrWN%.fa}_filteredPairedCov.txt
  
  java -jar ${progSplitSeqLowCov} -i ${mergedCorrWN%.fa}_filteredCov.txt -paired ${mergedCorrWN%.fa}_filteredPairedCov.txt -o ${mergedCorrWN%.fa}_filteredNotCov.txt -mCov 1 -fasta ${mergedCorrWN} -fastaOut ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  echo ${mergedCorrWN%.fa}_splitFiltered.fa >> $log
  java -jar ${progFastaStats} -i ${mergedCorrWN%.fa}_splitFiltered.fa -min 200 >> $log
  

# optional run RaGOO
########################################################## 

# 7. run AlignGraph for reference guided scallfolding 
########################################################## 
  echo "run AlignGraph..." 
  echo "run AlignGraph..." >> $log
  
  AlignGraphFolder=${workPath}/AlignGraph
  mkdir $AlignGraphFolder
  cd $AlignGraphFolder
  #convet fastq to fasta
  ${progseqtk} seq -a ${read1TrimPair[i]} > reads1.fa
  ${progseqtk} seq -a ${read2TrimPair[i]} > reads2.fa
  #run AlignGraph
  ${progAlignGraph} \
  --read1 reads1.fa --read2 reads2.fa \
  --contig ${mergedCorrWN%.fa}_splitFiltered.fa \
  --genome ${ref} \
  --distanceLow ${insLow} --distanceHigh ${insHight} --fastMap \
  --extendedContig extendedContigs.fa --remainingContig remainingContigs.fa \
  --misassemblyRemoval

  cat extendedContigs.fa remainingContigs.fa > AlignGraph_scallfolding.fa

# 8. run pilon 
########################################################## 
  echo "run pilon..." 
  echo "run pilon..." >> $log

  pilonFolder=${workPath}/pilon
  mkdir $pilonFolder
  cd $pilonFolder

  bwa index /lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/bin/RECORD/mod_PYGMY/results/edited_ref/edited_ref.fa

  bwa mem -t 20 /lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/bin/RECORD/mod_PYGMY/results/edited_ref/edited_ref.fa /lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/PYGMY/Trim_het/lib1_R1_trimPair.fastq /lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/PYGMY/Trim_het/lib1_R2_trimPair.fastq > reads.sam

  #only keep mapped reads
  samtools view -b -F 4 reads.sam -@ ${NThreads} > reads.bam

  samtools sort reads.bam -@ ${NThreads} > reads_sorted.all.bam

  samtools index reads_sorted.all.bam

  #######subset chromosomes
  CHR=`cat CHR.list`

  for i in $CHR;do

  samtools view -b -@ ${NThreads} reads_sorted.all.bam ${i} > ${i}.bam

  samtools index ${i}.bam

  java -jar -Xmx200G /lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/bin/pilon-1.23.jar --genome /lustre/nobackup/WUR/ABGC/liu194/analysis/SUS_PAN/Annotation/PYGMY/Split_chr/${i}.fa --frags ${i}.bam --output corrected_${i}.fasta --outdir ./ --changes --fix all --threads ${NThreads}

  done

###################################
  echo "FINISH!!!!!! ENJOY YOUR NEW GENOME ^_^" 
  echo "-. .-.   .-. .-.   .-. .-.   ."
  echo "  \   \ /   \   \ /   \   \ /"
  echo " / \   \   / \   \   / \   \"
  echo "~   \`-~ \`-\`   \`-~ \`-\`   \`-~ \`-"
  echo "FINISH!!!!!!  ^_^" >> $log
