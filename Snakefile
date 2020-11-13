#!/usr/bin/env python3
import pathlib2
import os
import pandas

#############
# FUNCTIONS #
#############

def resolve_path(x):
    return(str(pathlib2.Path(x).resolve(strict=False)))

def find_read_files(read_dir):
#Make list of files
    path_generator = os.walk(read_dir, followlinks = True)
    my_files = list((dirpath, filenames)
        for (dirpath, dirname, filenames)
        in path_generator)
#Make new dictionary & populate with files (flowcell = key)
    my_fastq_files = {}
    for dirpath, filenames in my_files:
        for filename in filenames:
            if filename.endswith('.fastq.gz'):
                my_flowcell = pathlib2.Path(dirpath).name
                my_fastq = str(pathlib2.Path(dirpath,filename))
                if my_flowcell in my_fastq_files:
                    my_fastq_files[my_flowcell].append(my_fastq)
                else:
                    my_fastq_files[my_flowcell]= []
                    my_fastq_files[my_flowcell].append(my_fastq)
    return(my_fastq_files)

def sample_name_to_fastq(wildcards):
    sample_row = sample_key[sample_key['Sample_name'] == wildcards.sample]
    sample_id = sample_row.iloc[-1]['OGF_sample_ID']
    sample_flowcell = sample_row.iloc[-1]['Flow_cell']
    sample_all_fastq = [x for x in all_fastq[sample_flowcell]
                        if '-{}-'.format(sample_id) in x]
    sample_r1 = sorted(list(x for x in sample_all_fastq
                            if '_R1_' in os.path.basename(x)))
    sample_r2 = sorted(list(x for x in sample_all_fastq
                            if '_R2_' in os.path.basename(x)))
    return({'r1': sample_r1, 'r2': sample_r2})

###########
# GLOBALS #
###########

read_dir = 'data/reads'

sample_key_file = 'data/sample_key.csv'

bbduk_adapters = '/adapters.fa'

star_reference_folder = 'output/star/star_reference'

#containers
bbduk_container = 'shub://TomHarrop/singularity-containers:bbmap_38.00'
salmon_container = 'shub://TomHarrop/singularity-containers:salmon_0.11.1'
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.7beta'

#########
# SETUP #
#########
# generate name to filename dictionary
all_fastq = find_read_files(read_dir)

sample_key = pandas.read_csv(sample_key_file)

all_samples = sorted(set(sample_key['Sample_name']))

#########
# RULES #
#########

rule target:
    input:
     #expand('output/asw_salmon/{sample}_quant/quant.sf', sample = all_samples),
     #'output/deseq2/ruakura/unann/unann_degs_blastx.outfmt6',
     #'output/deseq2/viral_expressed/no_annot/blastx.outfmt6',
     #'output/deseq2/ruakura/unann/interpro_degs.fasta.tsv',
     #'output/deseq2/viral_expressed/no_annot/interpro_degs.fasta.tsv',
     expand('output/asw_mh_concat_salmon/{sample}_quant/quant.sf', sample = all_samples),
     #expand('output/kraken/reports/kraken_{sample}_report.txt', sample=all_samples),
     'output/fastqc'

rule kraken:
    input:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz',
        db = 'data/20180917-krakendb'
    output:
        out = 'output/kraken/kraken_{sample}_out.txt',
        report = 'output/kraken/reports/kraken_{sample}_report.txt'
    log:
        'output/logs/kraken/kraken_{sample}.log'
    threads:
        20
    singularity:
        kraken_container
    shell:
        'kraken2 '
        '--threads {threads} '
        '--db {input.db} '
        '--paired '
        '--output {output.out} '
        '--report {output.report} '
        '--use-names '
        '{input.r1} {input.r2} '
        '&> {log}'

rule asw_mh_concat_salmon_quant:
    input:
        index_output = 'output/asw_mh_concat_salmon/transcripts_index/hash.bin',
        left = 'output/bbduk_trim/{sample}_r1.fq.gz',
        right = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/asw_mh_concat_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/asw_mh_concat_salmon/transcripts_index',
        outdir = 'output/asw_mh_concat_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon/asw_mh_concat_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.left} '
        '-2 {input.right} '
        '-o {params.outdir} '
        '--writeUnmappedNames '
        '-p {threads} '
        '&> {log}'

rule asw_mh_concat_salmon_index:
    input:
        transcriptome_length_filtered = 'data/asw_mh_transcriptome/asw_mh_isoforms_by_length.fasta'
    output:
        'output/asw_mh_concat_salmon/transcripts_index/hash.bin'
    params:
        outdir = 'output/asw_mh_concat_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/asw_mh_concat_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule interproscan_ru_behaviour_unann_degs:
	input:
		interpro_degs = 'output/deseq2/ruakura/unann/interpro_degs.fasta'
	output:
		interpro_tsv = 'output/deseq2/ruakura/unann/interpro_degs.fasta.tsv'
	params:
		outdir = 'output/deseq2/ruakura/unann'
	threads:
		20
	shell:
		'bin/interproscan-5.31-70.0/interproscan.sh '
		'--input {input.interpro_degs} '
		'--seqtype n '
		'--output-dir {params.outdir} '
		'--cpu {threads} '
		'--goterms'

rule filter_ru_behaviour_degs_for_interproscan:
    input:
        deg_ids = 'output/deseq2/ruakura/unann/interproscan_ids.txt',
        transcriptome_length_filtered = 'data/asw_transcriptome/isoforms_by_length.fasta'
    output:
       	interpro_degs = 'output/deseq2/ruakura/unann/interpro_degs.fasta'
    singularity:
        bbduk_container
    log:
        'filter_ru_behaviour_degs_for_interproscan'
    shell:
        'filterbyname.sh '
        'in={input.transcriptome_length_filtered} '
        'include=t '
        'substring=t '
        'names={input.deg_ids} '
        'out={output.interpro_degs} '
        '2> {log}'

rule interproscan_viral_unann_degs:
	input:
		interpro_degs = 'output/deseq2/viral_expressed/no_annot/interpro_degs.fasta'
	output:
		interpro_tsv = 'output/deseq2/viral_expressed/no_annot/interpro_degs.fasta.tsv'
	params:
		outdir = 'output/deseq2/viral_expressed/no_annot'
	threads:
		20
	shell:
		'bin/interproscan-5.31-70.0/interproscan.sh '
		'--input {input.interpro_degs} '
		'--seqtype n '
		'--output-dir {params.outdir} '
		'--cpu {threads} '
		'--goterms'

rule filter_viral_degs_for_interproscan:
    input:
        deg_ids = 'output/deseq2/viral_expressed/interproscan/interproscan_ids.txt',
        transcriptome_length_filtered = 'data/asw_transcriptome/isoforms_by_length.fasta'
    output:
       	interpro_degs = 'output/deseq2/viral_expressed/no_annot/interpro_degs.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_degs_for_interproscan.log'
    shell:
        'filterbyname.sh '
        'in={input.transcriptome_length_filtered} '
        'include=t '
        'substring=t '
        'names={input.deg_ids} '
        'out={output.interpro_degs} '
        '2> {log}'

rule blast_viral_unann_degs:
    input:
        unann_degs = 'output/deseq2/viral_expressed/no_annot/degs_no_annot.fasta'
    output:
        blastx_res = 'output/deseq2/viral_expressed/no_annot/blastx.outfmt6'
    params:
        blastdb = 'bin/blastdb/nr/nr'
    threads:
        40
    log:
        'output/logs/blast_viral_unann_degs.log'
    shell:
        'blastx '
        '-query {input.unann_degs} '
        '-db {params.blastdb} '
        '-num_threads {threads} '
        '-outfmt "6 std salltitles" > {output.blastx_res} '

rule filter_viral_degs_no_annot:
    input:
        deg_ids = 'output/deseq2/viral_expressed/no_annot/degs_no_annot.txt',
        transcriptome_length_filtered = 'data/asw_transcriptome/isoforms_by_length.fasta'
    output:
        unann_degs = 'output/deseq2/viral_expressed/no_annot/degs_no_annot.fasta'
    singularity:
        bbduk_container
    log:
        'output/logs/filter_viral_degs_no_annot.log'
    shell:
        'filterbyname.sh '
        'in={input.transcriptome_length_filtered} '
        'include=t '
        'substring=t '
        'names={input.deg_ids} '
        'out={output.unann_degs} '
        '2> {log}'

rule blast_ru_unann_genes:
    input:
        ru_unann_genes = 'output/deseq2/ruakura/unann/unann_degs.fasta'
    output:
        blast_res = 'output/deseq2/ruakura/unann/unann_degs_blastx.outfmt6'
    params:
        blastdb = 'bin/blastdb/nr/nr'
    threads:
        30
    log:
        'output/logs/blast_ru_unann_genes.log'
    shell:
        'blastx '
        '-query {input.ru_unann_genes} '
        '-db {params.blastdb} '
        '-num_threads {threads} '
        '-outfmt "6 std salltitles" > {output.blast_res}'

rule filter_ru_unann_genes:
    input:
        ru_unann_gene_ids = 'output/deseq2/ruakura/degs_with_no_annot.txt',
        asw_transcriptome = 'data/asw_transcriptome/isoforms_by_length.fasta'
    output:
        ru_unann_genes = 'output/deseq2/ruakura/unann/unann_degs.fasta'
    threads:
        50
    singularity:
        bbduk_container
    log:
        'output/logs/filter_ru_unann_genes.log'
    shell:
        'filterbyname.sh '
        'in={input.asw_transcriptome} '
        'include=t '
        'substring=t '
        'names={input.ru_unann_gene_ids} '
        'out={output.ru_unann_genes} '
        '&> {log}'

rule asw_salmon_quant:
    input:
        index_output = 'output/asw_salmon/transcripts_index/hash.bin',
        trimmed_r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        trimmed_r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    output:
        quant = 'output/asw_salmon/{sample}_quant/quant.sf'
    params:
        index_outdir = 'output/asw_salmon/transcripts_index',
        outdir = 'output/asw_salmon/{sample}_quant'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/salmon/asw_salmon_quant_{sample}.log'
    shell:
        'salmon quant '
        '-i {params.index_outdir} '
        '-l ISR '
        '-1 {input.trimmed_r1} '
        '-2 {input.trimmed_r2} '
        '-o {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule asw_salmon_index:
    input:
        transcriptome_length_filtered = 'data/asw_transcriptome/isoforms_by_length.fasta'
    output:
        'output/asw_salmon/transcripts_index/hash.bin'
    params:
        outdir = 'output/asw_salmon/transcripts_index'
    threads:
        20
    singularity:
        salmon_container
    log:
        'output/logs/asw_salmon_index.log'
    shell:
        'salmon index '
        '-t {input.transcriptome_length_filtered} '
        '-i {params.outdir} '
        '-p {threads} '
        '&> {log}'

rule fastqc:
    input:
        expand('output/bbduk_trim/{sample}_r{n}.fq.gz',
            sample=all_samples, n=[1,2])
    output:
        directory('output/fastqc')
    shell:
        'mkdir -p {output} ; '
        'fastqc --outdir {output} {input}'

rule bbduk_trim:
    input:
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    output:
        r1 = 'output/bbduk_trim/{sample}_r1.fq.gz',
        r2 = 'output/bbduk_trim/{sample}_r2.fq.gz'
    params:
        adapters = bbduk_adapters
    log:
        'output/logs/bbduk_trim/{sample}.log'
    threads:
        20
    singularity:
        bbduk_container
    shell:
        'bbduk.sh '
        'in={input.r1} '
        'in2={input.r2} '
        'out={output.r1} '
        'out2={output.r2} '
        'ref={params.adapters} '
        'ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=15 '
        '&> {log}'

##samples over 3 lanes this time - does this still work to concat all reads - think so?

rule cat_reads:
    input:
        unpack(sample_name_to_fastq)
    output: 
        r1 = temp('output/joined/{sample}_r1.fq.gz'),
        r2 = temp('output/joined/{sample}_r2.fq.gz')
    threads:
        1
    shell:
        'cat {input.r1} > {output.r1} & '
        'cat {input.r2} > {output.r2} & '
        'wait'