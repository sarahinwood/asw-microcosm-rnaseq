#!/usr/bin/env python3
import pathlib2
import os
import pandas
import peppy

#############
# FUNCTIONS #
#############

def get_reads(wildcards):
    input_keys = ['l6r1', 'l7r1', 'l8r1', 'l6r2', 'l7r2', 'l8r2']
    my_pep = pep.get_sample(wildcards.sample).to_dict()
    return {k: my_pep[k] for k in input_keys}

###########
# GLOBALS #
###########

##this parses the config & sample key files into an object named pep
pepfile: 'data/config.yaml'
##can now use this to generate list of all samples
all_samples = pep.sample_table['sample_name']
#make sample_table.csv in libre - excel format csv doesn't work
##to test it is working go into python interpreter inside venv
#import peppy
#proj = peppy.Project('path/to/config.yaml')
#proj.get_sample('sample_name')



bbduk_adapters = '/adapters.fa'
star_reference_folder = 'output/star/star_reference'

#containers
bbduk_container = 'shub://TomHarrop/seq-utils:bbmap_38.76'
salmon_container = 'docker://combinelab/salmon:latest'
kraken_container = 'shub://TomHarrop/singularity-containers:kraken_2.0.7beta'

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
     #'output/fastqc'

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
        index_output = 'output/asw_mh_concat_salmon/transcripts_index/refseq.bin',
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
        'output/asw_mh_concat_salmon/transcripts_index/refseq.bin'
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
        index_output = 'output/asw_salmon/transcripts_index/refseq.bin',
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
        'output/asw_salmon/transcripts_index/refseq.bin'
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

rule join_reads:
    input:
        unpack(get_reads)
    output: 
        r1 = 'output/joined/{sample}_r1.fq.gz',
        r2 = 'output/joined/{sample}_r2.fq.gz'
    threads:
        1
    shell:
        'cat {input.l6r1} {input.l7r1} {input.l8r1} > {output.r1} & '
        'cat {input.l6r2} {input.l7r2} {input.l8r2} > {output.r2} & '
        'wait'