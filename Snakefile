# Snakefile 
# Prerequisites (on Gadi): source /g/data/pq84/malaria/snakemake_pipeline/snakemake_setup.sh


################################################ Define paths using config file #########################################################################


configfile: "config/config.yaml"
fasta_path_pk=config['fasta_path_pk']
fasta_path_hg=config['fasta_path_hg']
bed_file_path=config['bed_file_path']
gff_path_pk=config['gff_path_pk']
gff_path_hg=config['gff_path_hg']
star_path=config['star_path']
source_dir=config['source_dir']
adapters=config['adapters']
rna_bloom=config['rna_bloom']
busco_db=config['busco_db']
busco_offline=config['busco_offline']

################################################ Set up: define input files, wildcards, etc. #########################################################################

## Defining the samples to be used for the {sample} wildcard
SAMPLES, = glob_wildcards(f'{source_dir}/{{sample}}_1.fq.gz') 

################################################## Define final files #######################################################################

rule all:
    input:
        "output/assembly/counts/abundance.tsv",
        "output/assembly/merged/merged-derep-assemblies.fa.idx",
        "output/assembly/merged/quality/missing_busco_list.tsv",
        "output/assembly/merged/merged-derep-assemblies.fa",
        "output/assembly/merged/merged-assemblies.fa.transdecoder.bed",
        "output/assembly/merged/merged-assemblies.fa",
        "output/trinity_assembly/trans_abundance/trinity.isoform.counts.matrix",
        "output/trinity_assembly/trans_abundance/trinity.gene.counts.matrix"

###################################################### Rules ##############################################################################

# Define local rules - not run with scheduler
localrules: all

# TRINITY ###########################################################################################################  

# Check quality of Trinity assembly - BUSCO
busco_trinity_pfx="output/trinity_assembly/quality/"

rule busco_quality:
    input:
        "output/trinity_assembly/Trinity.fasta"
    output:
        *[busco_trinity_pfx + filename for filename in [
            'run_plasmodium_odb10/full_table.tsv',
            'run_plasmodium_odb10/short_summary.json',
            'run_plasmodium_odb10/short_summary.txt',
            'missing_busco_list.tsv'
            ]]
    params:
        busco_trinity_pfx=f"{busco_trinity_pfx}",
        busco_db=busco_db,
        busco_offline=busco_offline
    shell:
        """
        busco -i {input} -l {params.busco_db} -m transcriptome -o {params.busco_trinity_pfx} -c 28 -f --offline --download_path {params.busco_offline}
        """

# Get isoform and gene abundances from trinity
trin_abun_pfx="output/trinity_assembly/trans_abundance/"

rule trinity_isoform_matrix:
    input:
        isoforms=expand(f"{trin_abun_pfx}{{sample}}/RSEM.isoforms.results", sample = SAMPLES),
        genes=expand(f"{trin_abun_pfx}{{sample}}/RSEM.genes.results", sample = SAMPLES),
        trinity_gene_map="output/trinity_assembly/Trinity.fasta.gene_trans_map"
    output:
        *[trin_abun_pfx + filename for filename in [
            'trinity.isoform.counts.matrix',
            'trinity.gene.counts.matrix'
            ]]
    params:
        trin_abun_pfx=f"{trin_abun_pfx}/trinity",
        input_list=lambda w: " ".join(expand(f"{trin_abun_pfx}{{sample}}/RSEM.isoforms.results", sample = SAMPLES))
    shell:
        """
        $TRINITY_HOME/util/abundance_estimates_to_matrix.pl \
        --est_method RSEM  \
        --cross_sample_norm TMM \
        --gene_trans_map {input.trinity_gene_map} \
        --out_prefix {params.trin_abun_pfx} \
        --name_sample_by_basedir \
        {params.input_list}
        """

rule trinity_abundances:
    input:
        trinity_fasta="output/trinity_assembly/Trinity.fasta",
        trinity_gene_map="output/trinity_assembly/Trinity.fasta.gene_trans_map",
        f_in="output/mapped/human/fixed_{sample}_Unmapped.out.mate1",
        r_in="output/mapped/human/fixed_{sample}_Unmapped.out.mate2",
        index="output/trinity_assembly/Trinity.fasta.bowtie2.rev.2.bt2"
    output:
        *[trin_abun_pfx + filename for filename in [
            '{sample}/RSEM.isoforms.results',
            '{sample}/RSEM.genes.results'
            ]]
    params:
        trin_abun_pfx=f"{trin_abun_pfx}{{sample}}"
    shell:
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl \
        --transcripts {input.trinity_fasta} \
        --seqType fq \
        --left {input.f_in} --right {input.r_in} \
        --est_method RSEM \
        --aln_method bowtie2 \
        --gene_trans_map {input.trinity_gene_map} \
        --output_dir {params.trin_abun_pfx} \
        --SS_lib_type RF \
        --thread_count 28
        """

trin_fast_pfx="output/trinity_assembly/"

rule trinity_abundances_prep:
    input:
        f"{trin_fast_pfx}Trinity.fasta"
    output:
        *[trin_fast_pfx + filename for filename in [
            'Trinity.fasta.bowtie2.1.bt2',
            'Trinity.fasta.bowtie2.2.bt2',
            'Trinity.fasta.bowtie2.3.bt2',
            'Trinity.fasta.bowtie2.4.bt2'
            'Trinity.fasta.bowtie2.rev.1.bt2',
            'Trinity.fasta.bowtie2.rev.2.bt2'
            ]]
    shell:
        """
        $TRINITY_HOME/util/align_and_estimate_abundance.pl \
        --transcripts {input} \
        --prep_reference \
        --est_method RSEM \
        --aln_method bowtie2 
        """

# Assemble with trinity
## reccomends 1GB of ram for every M reads - doubled it

rule trinity:
    input:
        f_in=expand("output/mapped/human/fixed_{sample}_Unmapped.out.mate1", sample = SAMPLES),
        r_in=expand("output/mapped/human/fixed_{sample}_Unmapped.out.mate2", sample = SAMPLES)
    output:
        "output/trinity_assembly/Trinity.fasta",
        "output/trinity_assembly/Trinity.fasta.gene_trans_map"
    params:
        left_list=lambda w: ",".join(expand("output/mapped/human/fixed_{sample}_Unmapped.out.mate1", sample = SAMPLES)),
        right_list=lambda w: ",".join(expand("output/mapped/human/fixed_{sample}_Unmapped.out.mate2", sample = SAMPLES)),
        out_pfx="output/trinity_assembly/"
    shell:
        """
        Trinity \
        --seqType fq \
        --max_memory 128G \
        --CPU 28 \
        --left {params.left_list} --right {params.right_list} \
        --SS_lib_type RF \
        --full_cleanup \
        --output {params.out_pfx}
        """

# Amend files for trinity

rule fix_header_for_trinity:
    input:
        f_in="output/mapped/human/{sample}_Unmapped.out.mate1",
        r_in="output/mapped/human/{sample}_Unmapped.out.mate2"
    output:
        f_out=temp("output/mapped/human/fixed_{sample}_Unmapped.out.mate1"),
        r_out=temp("output/mapped/human/fixed_{sample}_Unmapped.out.mate2")
    run:
        with open(input.f_in, "rt") as file_in:
            with open(output.f_out, "wt") as file_out:
                for line in file_in:
                    file_out.write(line.replace('0:N', '1:N'))

        with open(input.r_in, "rt") as file_in:
            with open(output.r_out, "wt") as file_out:
                for line in file_in:
                    file_out.write(line.replace('1:N', '2:N'))

# RNA BLOOM ########################################################################################################### 

## Quantify
kallisto_pfx="output/assembly/counts/"

def zip_my_files(w, input):
    a = input.sample_mates
    return sum([list(t) for t in list(zip(a[:len(a)//2 + 1], a[len(a)//2:]))], [])

rule kallisto_quant:
    localrule: True
    input:
        index="output/assembly/merged/merged-derep-assemblies.fa.idx",
        sample_mates=expand(["output/mapped/human/{sample}_Unmapped.out.mate1", "output/mapped/human/{sample}_Unmapped.out.mate2"], sample= SAMPLES),
    output:
        *[kallisto_pfx + filename for filename in [
            'abundance.h5',
            'abundance.tsv',
            'run_info.json'
            ]]
    params:
        kallisto_pfx=kallisto_pfx,
        sample_list = zip_my_files,
    shell:
        """
        kallisto quant -o {params.kallisto_pfx} -t 10 -i {input.index} {params.sample_list}
        """

## Index
rule kallisto_index_transcriptome:
    input:
        "output/assembly/merged/merged-derep-assemblies.fa"
    output:
        "output/assembly/merged/merged-derep-assemblies.fa.idx"
    shell:
        """
        kallisto index -i {output} -k 31 -t 5 {input}
        """

# Check quality of assembly - BUSCO
busco_pfx="output/assembly/merged/quality/"

rule busco_quality:
    input:
        "output/assembly/merged/merged-derep-assemblies.fa"
    output:
        *[busco_pfx + filename for filename in [
            'run_plasmodium_odb10/full_table.tsv',
            'run_plasmodium_odb10/short_summary.json',
            'run_plasmodium_odb10/short_summary.txt',
            'missing_busco_list.tsv'
            ]]
    params:
        busco_pfx=f"{busco_pfx}",
        busco_db=busco_db,
        busco_offline=busco_offline
    shell:
        """
        busco -i {input} -l {params.busco_db} -m transcriptome -o {params.busco_pfx} -c 28 -f --offline --download_path {params.busco_offline}
        """

# Remove duplicates
rule cd_hit_clean:
    input:
        "output/assembly/merged/merged-assemblies.fa"
    output:
        "output/assembly/merged/merged-derep-assemblies.fa"
    shell:
        """
        cd-hit-est -i {input} -o {output} -c 0.95 -n 10 -p 1 -g 1 -M 600000 -T 5 -d 40 
        """

# Annotate open reeding frames
merge_pfx="output/assembly/merged/"

rule transdecoder_coding_regions:
    input:
        "output/assembly/merged/merged-assemblies.fa"
    output:
        *[merge_pfx + filename for filename in [
            'merged-assemblies.fa.transdecoder.pep',
            'merged-assemblies.fa.transdecoder.cds',
            'merged-assemblies.fa.transdecoder.gff3',
            'merged-assemblies.fa.transdecoder.bed'
            ]]
    params:
        merge_pfx
    shell:
        """
        TransDecoder.LongOrfs -t {input} --output_dir {params} -m 150
        TransDecoder.Predict -t {input} --single_best_only --output_dir {params} 
        """

# Merge RNA bloom assemblies with transabyss

rule transabyss_merge_assemblies:
    input:
        expand("output/assembly/{sample}/rnabloom.transcripts.fa", sample = SAMPLES)
    output:
        "output/assembly/merged/merged-assemblies.fa"
    params:
        fasta=lambda w: " ".join(expand("output/assembly/{sample}/rnabloom.transcripts.fa", sample = SAMPLES)),
        pfx=lambda w: ". ".join(expand("{sample}", sample = SAMPLES))
    shell:
        """
        transabyss-merge \
        --mink 32 --maxk 32 \
        --out {output} \
        --threads 5 \
        {params.fasta} \
        --prefix {params.pfx}
        """
        
# Assemble with RNA bloom

rule rna_bloom:
    input:
        rna_bloom=rna_bloom,
        f_in="output/mapped/human/{sample}_Unmapped.out.mate1",
        r_in="output/mapped/human/{sample}_Unmapped.out.mate2"
    output:
        "output/assembly/{sample}/rnabloom.transcripts.fa",
        "output/assembly/{sample}/rnabloom.transcripts.short.fa",
        "output/assembly/{sample}/rnabloom.transcripts.nr.fa"
    params:
        "output/assembly/{sample}"
    shell:
        """
        java -jar {input.rna_bloom} \
        -left {input.f_in} -right {input.r_in} \
        -revcomp-right \
        -t 21 \
        -outdir {params}
        """

# QC ###########################################################################################################  

# Map to Pk with STAR using unmapped reads post mapping to hg

## Max intron size of Pf is ~1000 & min is 5- specified in input

pfx_pk="output/mapped/"

rule star_aligner_pk:
    input:
        "output/genome_index/Genome",
        "output/genome_index/SA",
        "output/genome_index/SAindex",
        star=star_path,
        f_in="output/mapped/human/{sample}_Unmapped.out.mate1",
        r_in="output/mapped/human/{sample}_Unmapped.out.mate2"
    output:
        *[pfx_pk + filename for filename in [
            '{sample}_Aligned.sortedByCoord.out.bam',
            '{sample}_Log.final.out',
            '{sample}_Log.out',
            '{sample}_Log.progress.out',
            '{sample}_SJ.out.tab',
            '{sample}_ReadsPerGene.out.tab'
            ]]
    params:
        pfx_pk=f"{pfx_pk}{{sample}}_",
        genome="output/genome_index"
    shell:
        """
        {input.star} \
        --runThreadN 5 \
        --genomeDir {params.genome} \
        --readFilesIn {input.f_in} {input.r_in} \
        --outSAMtype BAM SortedByCoordinate \
	    --outFileNamePrefix {params.pfx_pk} \
	    --outSAMattributes Standard \
        --alignIntronMin 5 \
        --limitBAMsortRAM 1847354302 \
        --alignIntronMax 1000 \
        --quantMode GeneCounts
        """

# Map to hg with STAR

## Max intron size of Pf is ~1000 & min is 5- specified in input

pfx_hg="output/mapped/human/"

rule star_aligner_hg:
    input:
        "output/genome_index/human/Genome",
        "output/genome_index/human/SA",
        "output/genome_index/human/SAindex",
        star=star_path,
        f_in="output/trimmed/{sample}_trimmed_1.fq.gz",
        r_in="output/trimmed/{sample}_trimmed_2.fq.gz"
    output:
        *[pfx_hg + filename for filename in [
            '{sample}_Aligned.sortedByCoord.out.bam',
            '{sample}_Log.final.out',
            '{sample}_Log.out',
            '{sample}_Log.progress.out',
            '{sample}_SJ.out.tab',
            '{sample}_ReadsPerGene.out.tab',
            '{sample}_Unmapped.out.mate1',
            '{sample}_Unmapped.out.mate2'
            ]]
    params:
        pfx_hg=f"{pfx_hg}{{sample}}_",
        genome="output/genome_index/human"
    shell:
        """
        {input.star} \
        --runThreadN 5 \
        --genomeDir {params.genome} \
        --readFilesIn {input.f_in} {input.r_in} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
	    --outFileNamePrefix {params.pfx_hg} \
	    --outSAMattributes Standard \
        --alignIntronMin 5 \
        --limitBAMsortRAM 1847354302 \
        --alignIntronMax 1000 \
        --quantMode GeneCounts \
        --outReadsUnmapped Fastx 
        """

# Index genomes for mapping STAR

## Pk

pfx_pk_index="output/genome_index/"

rule star_index_pk:
    input:
        fasta_pk=fasta_path_pk,
        gff_pk=gff_path_pk,
        star=star_path
    output:
        *[pfx_pk_index + filename for filename in [
            'Genome','SA','SAindex','chrLength.txt','chrName.txt',
            'chrStart.txt','exonGeTrInfo.tab','exonInfo.tab',
            'geneInfo.tab','genomeParameters.txt','jdbInfo.txt',
            'sjdbList.fromGTF.out.tab','sjdbList.out.tab',
            'transcriptInfo.tab'
            ]]
    params:
    params:
        "output/genome_index"
    shell:
        """
        {input.star} \
        --runThreadN 5 \
        --runMode genomeGenerate \
        --genomeSAindexNbases 11 \
        --genomeDir {params} \
        --genomeFastaFiles {input.fasta} \
        --sjdbGTFfile {input.gff} \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon exon \
        --sjdbOverhang 149 \
        --outFilterMultimapNmax 20
        """

# Human

pfx_hg_index="output/genome_index/human/"

rule star_index_hg:
    input:
        fasta_hg=fasta_path_hg,
        gff_hg=gff_path_hg,
        star=star_path
    output:
        *[pfx_hg_index + filename for filename in [
            'Genome','SA','SAindex','chrLength.txt','chrName.txt',
            'chrStart.txt','exonGeTrInfo.tab','exonInfo.tab',
            'geneInfo.tab','genomeParameters.txt','jdbInfo.txt',
            'sjdbList.fromGTF.out.tab','sjdbList.out.tab',
            'transcriptInfo.tab'
            ]]
    params:
        "output/genome_index/human/"
    shell:
        """
        {input.star} \
        --runThreadN 5 \
        --runMode genomeGenerate \
        --genomeSAindexNbases 11 \
        --genomeDir {params} \
        --genomeFastaFiles {input.fasta_hg} \
        --sjdbGTFfile {input.gff_hg} \
        --sjdbGTFtagExonParentTranscript Parent \
        --sjdbGTFfeatureExon exon \
        --sjdbOverhang 149 \
        --outFilterMultimapNmax 20
        """


# Remove adapaters

rule bbduk:
    input:
        adapters=adapters,
        f_in=f"{source_dir}/{{sample}}_1.fq.gz",
        r_in=f"{source_dir}/{{sample}}_2.fq.gz"
    output:
        f_out="output/trimmed/{sample}_trimmed_1.fq.gz",
        r_out="output/trimmed/{sample}_trimmed_2.fq.gz",
        stats="output/trimmed/{sample}.stats"
    shell:
        """
        bbduk.sh \
        in1={input.f_in} in2={input.r_in} \
        out1={output.f_out} out2={output.r_out} \
        ktrim=r \
        k=12 \
        mink=6 \
        hdist=1 \
        ref={input.adapters} \
        tpe \
        tbo \
        stats={output.stats}
        """

# QC fastq files

rule qc:
    input:
        source_dir
    output:
        fast="output/fastqc/",
        multi="output/multiqc/"
    params:
        in_dir=f"{source_dir}/*fq.gz",
    shell:
        """
        mkdir output/fastqc

        fastqc -t 5 -o {output.fast} {params.in_dir}
        """