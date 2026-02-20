// ### UMI a read preprocessing ###
// # adaptors to trim: AACCGCCAGGAGT
// # UMI 1-8 bases in read position
process UMI_extract {
	tag "UMI_extract on $name using $task.cpus CPUs and $task.memory memory"
	//label "s_cpu"
	//label "xxs_mem"
	//container "quay.io/biocontainers/mulled-v2-452184f0fcbe6b0806405adb8f0b3873a7bd70a8:9c5efc89686d718e262836409bbf8541c6c5a159-0"

	input:
	tuple val(name), val(sample), path(fwd), path(rev)

	output:
	tuple val(name), val(sample), path("${name}.umi*.fastq.gz")

	script:
	"""
	echo UMI_extract $name
	source /home/anna/miniconda3/bin/activate umi_tools
	umi_tools extract -I ${fwd} --bc-pattern=NNNNNNNN --read2-in=${rev} --stdout=${name}.umi.R1.fastq.gz --read2-out=${name}.umi.R2.fastq.gz || \
	umi_tools extract --ignore-read-pair-suffixes -I ${fwd} --bc-pattern=NNNNNNNN --read2-in=${rev} --stdout=${name}.umi.R1.fastq.gz --read2-out=${name}.umi.R2.fastq.gz
	"""
}

	// --ignore-read-pair-suffixes // use this for MGI

process TRIMMING {
	tag "TRIMMING on $name using $task.cpus CPUs and $task.memory memory"
	//label "s_cpu"
	//label "xxs_mem"
	
	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("*.fastq.gz")

	script:
	"""
	echo TRIMMING $name
	source /home/anna/miniconda3/bin/activate cutadapt
	cutadapt -g AACCGCCAGGAGT -m 50 -o ${name}.trimmed1.R1.fastq.gz -p ${name}.trimmed1.R2.fastq.gz $reads
	"""
}


process FIRST_ALIGN_BAM {
	tag "FIRST_ALIGN_BAM on $name using $task.cpus CPUs and $task.memory memory"
	//label "m_cpu"
	//label "l_mem"

	input:
	tuple val(name), val(sample), path(reads)

	output:
	tuple val(name), val(sample), path("${name}.bam")

	script:
	rg = "\"@RG\\tID:${name}\\tSM:${name}\\tLB:${name}\\tPL:ILLUMINA\""
	"""
	echo FIRST_ALIGN_BAM $name
	source /home/anna/miniconda3/bin/activate bwa
	bwa mem -R ${rg} -t $task.cpus ${params.refindex} $reads > ${name}.sam
	samtools view -Sb ${name}.sam -o ${name}.bam
	"""
}

process SORT_INDEX {
	tag "Sort index on $name using $task.cpus CPUs and $task.memory memory"
	//label "m_mem"
	//label "xs_cpu"

	input:
	tuple val(name), val(sample), path(bam)

	output:
	tuple val(name), val(sample), path("${name}.sorted.bam"), path("${name}.sorted.bai")


	script:
	"""
	echo SORT_INDEX $name
	source /home/anna/miniconda3/bin/activate samtools
	samtools sort $bam -o ${name}.sorted.bam
	samtools index ${name}.sorted.bam ${name}.sorted.bai
	"""
}

process FLAGSTAT {
	tag "FLAGSTAT on $name"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.flagstat.txt")

	script:
	"""
	echo FLAGSTAT $name
	source /home/anna/miniconda3/bin/activate bwa-samtools
	samtools flagstat $bam > ${name}.flagstat.txt
	"""
}


process DEDUP {
	tag "DEDUP on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/mapped/", mode:'copy'
	//label "s_cpu"
	//label "m_mem"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.deduped.bam"), path("${name}.deduped.bai")

	script:
	"""
	echo DEDUP $name
	source /home/anna/miniconda3/bin/activate umi_tools
	umi_tools dedup -I $bam --paired -S ${name}.deduped.bam
	samtools index ${name}.deduped.bam ${name}.deduped.bai
	"""
}

process MUTECT2 {
	tag "MUTECT2 on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/vcfs/", mode:'copy'
	//label "m_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam), path(bai)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.vcf")

	script:
	"""
	echo MUTECT2 $name
	source /home/anna/miniconda3/bin/activate gatk4
	gatk Mutect2 --reference ${params.ref}.fa --input ${bam} --tumor-sample $name --annotation StrandArtifact --min-base-quality-score 20 --intervals $params.ivl --output ${sample.name}.mutect.vcf
	"""
}
//source /home/anna/miniconda3/bin/activate gatk400

process FILTER_MUTECT {
	tag "FILTER_MUTECT on $name using $task.cpus CPUs and $task.memory memory"
	//label "s_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.filt.vcf")

	script:
	"""
	echo FILTER_MUTECT $name
	source /home/anna/miniconda3/bin/activate gatk4
	gatk FilterMutectCalls -V $vcf_input -O ${sample.name}.mutect.filt.vcf
	"""
}

process NORMALIZE_MUTECT {
	tag "NORMALIZE_MUTECT on $name using $task.cpus CPUs $task.memory"
	//label "xxs_mem"
	//label "s_cpu"
	//container "staphb/bcftools:1.10.2"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path ("${sample.name}.mutect.filt.norm.vcf")

	script:
	"""
	echo NORMALIZE_MUTECT $name
    source /home/anna/miniconda3/bin/activate bcftools117
	bcftools norm -m-both $vcf_input > ${sample.name}.mutect.filt.norm.vcf
	"""
}

process ANNOTATE_MUTECT {
	tag "ANNOTATE_MUTECT on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/vcfs/", mode:'copy'
	//label "s_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${sample.name}.mutect2.filt.norm.vep.vcf")

	script:
	"""
	echo ANNOTATE_MUTECT $name
	source /home/anna/miniconda3/bin/activate  vep
	vep -i $vcf_input --cache --cache_version 95 --dir_cache $params.vep \
	--fasta ${params.ref}.fa --merged --offline --vcf --everything -o ${sample.name}.mutect2.filt.norm.vep.vcf
	"""	
}

process FILTER_VCF {
	tag "FILTER_VCF on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/vcfs/", mode:'copy'
	//container "staphb/bcftools:1.10.2"
	//label "xxs_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(name), val(sample), path("${sample.name}.mutect2.filt.norm.vep.filt.vcf")

	script:
	"""
	echo FILTER_VCF $name
    source /home/anna/miniconda3/bin/activate bcftools117
	bcftools view -f 'PASS,clustered_events,multiallelic' $vcf_input > ${sample.name}.mutect2.filt.norm.vep.filt.vcf
	"""	
}

process VCF2CSV {
	tag "VCF2CSV on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/variants/", mode:'copy'
	//label "xs_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(vcf_input)
	
	output:
	tuple val(sample.run), path("${name}.csv")

	script:
	"""
	echo VCF2CSV $name
	source /home/anna/miniconda3/bin/activate vcf2csv
	python $params.vcf2csv simple --build GRCh37 -i $vcf_input -o ${name}.csv
	"""	
}

process MERGE_TABLES {
	tag "MERGE_TABLES on $run using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${run}/variants/", mode:'copy'
	//label "s_mem"
	//label "s_cpu"
	debug true

	input:
	tuple val(run), path(all_annotated_normed)
	
	output:
	path "${run}.merged_variant.table_NEW.tsv"

	script:
	"""
	echo MERGE_TABLES $run
	source /home/anna/miniconda3/bin/activate erko
	Rscript --vanilla $params.mergescript $run
	"""	
}

process FLT3 {
	tag "FLT3 on $name using $task.cpus CPUs $task.memory"
	publishDir "${params.outDirectory}/${sample.run}/FLT3/", mode:'copy'
	//label "m_mem"
	//label "s_cpu"
	errorStrategy { task.exitStatus in [143,137,104,134,139,247,null,'', '-'] ? 'retry' : 'ignore' } //this does not really work

	input:
	tuple val(name), val(sample), path(bam), path(bai)
	
	output:
	tuple val(name), val(sample), path("${name}.deduped_FLT3_ITD_summary.txt")

	script:
	"""
	echo FLT3 $name
	source /home/anna/miniconda3/bin/activate perl
	tar -C /tmp -xf $params.flt3tar
	perl /tmp/FLT3/FLT3_ITD_ext/FLT3_ITD_ext.pl --bam $bam --output ./ --ngstype amplicon --genome hg19 --fgbiojar $params.fgbio --picardjar $params.picard --refindex /tmp/FLT3/FLT3_bwaindex/FLT3_dna_e1415
	touch ${name}.deduped_FLT3_ITD_summary.txt
	ls -alh
	"""	
}


process BAMQC {
	tag "BAMQC on $name using $task.cpus CPUs and $task.memory memory"
	//label "s_mem"
	//label "s_cpu"
	//container 'quay.io/biocontainers/mulled-v2-b0664646864bfdb46c5343b1b2b93fc05adb4b77:39a005770a3e30fb6aa3bf424b57ddf52bae7ece-0'

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("*")

	script:
	"""
	echo BAMQC $name
    source /home/anna/miniconda3/bin/activate qc_picard_samtools
	samtools flagstat $bam > ${name}.flagstat
	samtools stats $bam > ${name}.samstats
	picard CollectHsMetrics I=$bam BAIT_INTERVALS=$params.ivl TARGET_INTERVALS=$params.ivl R=${params.ref}.fa O=${name}.hs_metrics
	picard CollectAlignmentSummaryMetrics I=$bam R=${params.ref}.fa O=${name}.aln_metrics
	"""
}

process COVERAGE {
	tag "COVERAGE on $name using $task.cpus CPUs and $task.memory memory"
	//label "xl_mem"
	//label "s_cpu"

	input:
	tuple val(name), val(sample), path(bam), path(bai)

	output:
	tuple val(name), val(sample), path("${name}.PBcov.cons.txt")

	script:
	"""
	echo COVERAGE $name
	source /home/anna/miniconda3/bin/activate bedtools
	bedtools coverage -abam $params.bed -b $bam -d > ${name}.PBcov.cons.txt
	"""
}

process COVERAGE_POSTPROCESS {
	tag "COVERAGE_POSTPROCESS on $name using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${sample.run}/Cov/", mode:'copy'
	//label "s_mem"
	//label "xxs_mem"

	input:
	tuple val(name), val(sample), path(txt)

	output:
	tuple val(name), val(sample), path("${name}.perexon_stat.txt")

	script:
	"""
	echo COVERAGE_POSTPROCESS $name
	source /home/anna/miniconda3/bin/activate erko
	Rscript --vanilla $params.covscript $txt
	ls -al
	"""
}

process MULTIQC {
	tag "MultiQC on all samples from $run using $task.cpus CPUs and $task.memory memory"
	publishDir "${params.outDirectory}/${run}/QC/", mode:'copy'
	//container 'ewels/multiqc:v1.18'
	//label "xs_mem"
	//label "s_cpu"

	input:
	tuple val(run), path(all_annotated_normed)

	output:
	path "*"

	script:
	"""
	echo MULTIQC $run
	source /home/anna/miniconda3/bin/activate multiqc18
	multiqc . -n MultiQC.html
	"""

}

//	/home/anna/miniconda3/bin/multiqc . -n MultiQC.html


workflow {
	rawFQs = Channel.fromPath("${params.homeDir}/samplesheet.csv")
    .splitCsv(header: true)
    .map { row ->
        
        //def fileR1 = file("${params.inputDirectory}/${row.name}_R1.fastq.gz", checkIfExists: true)
        //def fileR2 = file("${params.inputDirectory}/${row.name}_R2.fastq.gz", checkIfExists: true)
		def fileR1 = file("${params.inputDirectory}/${row.run}/${row.name}_R1.fastq.gz", checkIfExists: true)
		def fileR2 = file("${params.inputDirectory}/${row.run}/${row.name}_R2.fastq.gz", checkIfExists: true)


//        def fileSizeR1 = fileR1.size()
//        def fileSizeR2 = fileR2.size()
//		def size = MemoryUnit.of(fileSizeR1 + fileSizeR2).toMega() / 1000

	def meta = [name: row.name, run: row.run]
//		def meta = [name: row.name, run: row.run, size: size]
        [
            meta.name,
            meta,
            fileR1,
            fileR2,
		]
    }
    .view()

	umiFQs = UMI_extract(rawFQs)
	trimmedFQs = TRIMMING(umiFQs)
	firstBAM = FIRST_ALIGN_BAM(trimmedFQs)
	sortedBamBai = SORT_INDEX(firstBAM)
	dedupedBam = DEDUP(sortedBamBai)
	QCs = BAMQC(dedupedBam)
	ToMultiQC = QCs.map({return [it[1].run, it[2]]}).
		groupTuple()
		.map({return [it[0], it[1].flatten()]})//.view{"$it is groupTuple"}
	MULTIQC(ToMultiQC)

	PBcov = COVERAGE(dedupedBam)
	COVERAGE_POSTPROCESS(PBcov)

	rawVCFs = MUTECT2(dedupedBam)
	filteredVCFs = FILTER_MUTECT(rawVCFs)
	normedVCFs = NORMALIZE_MUTECT(filteredVCFs)
	annotatedVCFs = ANNOTATE_MUTECT(normedVCFs)
	filteredAnnotatedVCFs = FILTER_VCF(annotatedVCFs)
	CSVs = VCF2CSV(filteredAnnotatedVCFs)
	MERGE_TABLES(CSVs.groupTuple())

	FLT3(dedupedBam)	
}
