import glob
import os
import pandas as pd


"""
this file performs the analyses for the odontosyllis luciferase
"""

#   vvvv  LOOK HERE vvv
#
#  your path to the directory containing trimmomatic-0.35.jar and adapters/TruSeq3-PE-2.fa
#   goes here
trimjar = "/usr/local/bin/Trimmomatic-0.35/trimmomatic-0.35.jar"
adapterpath = "TruSeq3-PE-2-fastqc.fa"
configfile: "config.yaml"
maxthreads = 90
dockeruser = "dschultz"
nrpath = "/data/ncbi/db/nr"
#
#
#   ^^^^  LOOK HERE ^^^

def sample2rna_f(wildcards):
    return config["species"][wildcards.taxon]["rna_f"]

def sample2rna_r(wildcards):
    return config["species"][wildcards.taxon]["rna_r"]

rule all:
    input:
        "blast_results/orthology_report.txt",
        expand("figures/{sample}_transcript_figure.pdf", sample = config["samples"])

rule download:
    output:
        "reads/{taxon}_f.fastq.gz",
    params:
        link=sample2rna_f,
    threads:
        maxthreads
    shell:
        "wget {params.link} -O {output}"

rule download_rev:
    output:
        "reads/{taxon}_r.fastq.gz",
    params:
        link=sample2rna_r,
    threads:
        maxthreads
    shell:
        "wget {params.link} -O {output}"

rule transcript_figure:
    input:
        reference = "fasta/isoform_analysis/c9g1i6.fasta",
        pepbam = "alignments/{sample}_pepfrags_to_c9g1i6_bwa.sorted.bam",
        annotation = "annotations/DN46871_c9g1i6.gff",
        #rnaalign = "alignments/{sample}_RNA_alignment_to_c9g1i6.sorted.bam",
        rnaalign = "alignments/{sample}_RNA_alignment_to_c9g1i6_bwa.sorted.bam",
        longcdna = "alignments/{sample}_longcdna_to_c9g1i6.sorted.bam"
    output:
        figname = "figures/{sample}_transcript_figure.pdf"
    shell:
        """pauvre browser -r {input.reference} \
        -c DN46871_c9g1i6 --start 1 --stop 1231 \
        --path {output.figname} \
        -p "bam:{input.pepbam}:reads:narrow,remdups" "gff3:{input.annotation}" "bam:{input.longcdna}:reads" "bam:{input.rnaalign}:depth:c" """

rule trim_polychaete_RNAseq:
    """This should have some failsafe to make sure that the trimmed reads
    actually have something in them..."""
    input:
        f1 = "reads/{taxon}_f.fastq.gz",
        f2 = "reads/{taxon}_r.fastq.gz",
        trim_jar = trimjar,
        adapter_path = adapterpath
    output:
        f_paired =   temp("reads/trimmed/{taxon}_polyRNAseq_f.trimmed.fastq.gz"),
        r_paired =   temp("reads/trimmed/{taxon}_polyRNAseq_r.trimmed.fastq.gz"),
        f_unpaired = temp("reads/trimmed/{taxon}_1.trim.unpaired.fastq.gz"),
        r_unpaired = temp("reads/trimmed/{taxon}_2.trim.unpaired.fastq.gz")
        #f_paired =   expand("reads/trimmed/{taxon}_polyRNAseq_f.trimmed.fastq.gz", sample = config["species"]),
        #r_paired =   expand("reads/trimmed/{taxon}_polyRNAseq_r.trimmed.fastq.gz", sample = config["species"]),
        #f_unpaired = temp(expand("reads/trimmed/{taxon}_1.trim.unpaired.fastq.gz", sample = config["species"])),
        #r_unpaired = temp(expand("reads/trimmed/{taxon}_2.trim.unpaired.fastq.gz", sample = config["species"]))

    threads:
        maxthreads - 1
    shell:
        """java -jar {input.trim_jar} PE \
        -threads {threads} \
        -phred33 {input.f1} {input.f2} \
        {output.f_paired} \
        {output.f_unpaired} \
        {output.r_paired} \
        {output.r_unpaired} \
        ILLUMINACLIP:{input.adapter_path}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36"""

rule rename_fastqs:
    input:
        f_paired = "reads/trimmed/{taxon}_polyRNAseq_f.trimmed.fastq.gz",
        r_paired = "reads/trimmed/{taxon}_polyRNAseq_r.trimmed.fastq.gz",
    output:
        f_renamed = protected("reads/trimmed/{taxon}_polyRNAseq_f.trimmed.renamed.fastq.gz"),
        r_renamed = protected("reads/trimmed/{taxon}_polyRNAseq_r.trimmed.renamed.fastq.gz"),
    shell:
        """zcat {input.f_paired} | \
        awk '{{print (NR%4 == 1) ? "@1_" ++i "/1": $0}}' | \
        gzip -c > {output.f_renamed};

        zcat {input.r_paired} | \
        awk '{{print (NR%4 == 1) ? "@1_" ++i "/2": $0}}' | \
        gzip -c > {output.r_renamed}"""

rule assemble_polychaete_txome:
    input:
        f_paired = "reads/trimmed/{taxon}_polyRNAseq_f.trimmed.renamed.fastq.gz",
        r_paired = "reads/trimmed/{taxon}_polyRNAseq_r.trimmed.renamed.fastq.gz",
    output:
        outpath = temp("txomes/trinity_{taxon}"),
        assemblypath = protected("txomes/trinity_{taxon}.fasta")
    params:
        outpath = "txomes/trinity_{taxon}",
        outfasta = "txomes/trinity_{taxon}/Trinity.fasta",
        dockeruser = dockeruser
    threads:
        maxthreads
    shell:
        """docker run \
        -u $(id -u):$(id -g) --rm  \
        -v `pwd`:`pwd` \
        -v /etc/passwd:/etc/passwd \
        trinityrnaseq/trinityrnaseq Trinity \
        --seqType fq \
        --left `pwd`/{input.f_paired} \
        --right `pwd`/{input.r_paired} \
        --max_memory 500G \
        --CPU {threads} \
        --output `pwd`/{params.outpath};
        cp {params.outfasta} {output.assemblypath}"""

rule txome_transdecoder:
    input:
        assemblypath = "txomes/trinity_{taxon}.fasta"
    output:
        transdecoder_peps = protected("txomes/pep/{taxon}_longest_orfs.pep"),
        transdecoder_dir = "txomes/pep/trinity_{taxon}.fasta.transdecoder_dir",
        temp_td_dir = temp("trinity_{taxon}.fasta.transdecoder_dir.__checkpoints_longorfs")
    threads:
        1
    params:
        temp_transdecoder_peps = "trinity_{taxon}.fasta.transdecoder_dir/longest_orfs.pep",
        temp_transdecoder_dir = "trinity_{taxon}.fasta.transdecoder_dir"
    shell:
        """transdecoder -t {input.assemblypath}; \
        cp {params.temp_transdecoder_peps} {output.transdecoder_peps}; \
        mv {params.temp_transdecoder_dir} txomes/pep/"""

rule make_protein_blastdb:
    """this makes protein blastdbs for all of the transcriptomes"""
    input:
        transdecoder_peps = "txomes/pep/{taxon}_longest_orfs.pep"
    output:
        phr = "txomes/pep/{taxon}_pep.phr",
        pin = "txomes/pep/{taxon}_pep.pin",
        pog = "txomes/pep/{taxon}_pep.pog",
        psd = "txomes/pep/{taxon}_pep.psd",
        psi = "txomes/pep/{taxon}_pep.psi",
        psq = "txomes/pep/{taxon}_pep.psq"
    threads:
        1
    params:
        dbname = "{taxon}_pep",
        outname = "txomes/pep/{taxon}_pep"
    shell:
        """makeblastdb -in {input.transdecoder_peps} \
        -dbtype prot \
        -title {params.dbname} \
        -parse_seqids \
        -out {params.outname}"""

rule make_nucl_blastdb:
    """this makes nucl blastdbs for all of the transcriptomes"""
    input:
        fasta_txome = "txomes/trinity_{taxon}.fasta"
    output:
        nhr = "txomes/{taxon}_nucl.nhr",
        nin = "txomes/{taxon}_nucl.nin",
        nog = "txomes/{taxon}_nucl.nog",
        nsd = "txomes/{taxon}_nucl.nsd",
        nsi = "txomes/{taxon}_nucl.nsi",
        nsq = "txomes/{taxon}_nucl.nsq"
    threads:
        1
    params:
        dbname = "{taxon}_nucl",
        outname = "txomes/{taxon}_nucl"
    shell:
        """makeblastdb -in {input.fasta_txome} \
        -dbtype nucl \
        -title {params.dbname} \
        -parse_seqids \
        -out {params.outname}"""

rule tblastn_c9g1i2:
    """this rule uses blastp to match the c9g1i2 transcript against
    the transcriptomes of interest. it outputs blast results to /blast_results
    """
    input:
        query = "fasta/luc_proteins/c9g1i2.pep",
        nhr = "txomes/{taxon}_nucl.nhr",
        nin = "txomes/{taxon}_nucl.nin",
        nog = "txomes/{taxon}_nucl.nog",
        nsd = "txomes/{taxon}_nucl.nsd",
        nsi = "txomes/{taxon}_nucl.nsi",
        nsq = "txomes/{taxon}_nucl.nsq"
    output:
        results = "blast_results/tblastn/c9g1i2_to_{taxon}_tblastn.txt"
    threads:
        maxthreads
    params:
        db = "txomes/{taxon}_nucl"
    shell:
        """tblastn -db {params.db} \
        -query {input.query} \
        -outfmt 6 \
        -db_gencode 1 \
        -max_target_seqs 1 \
        -num_threads {threads} > {output.results}"""

rule get_top_hits_tblastn:
    """ gets the nucl sequence of the top hit out of the target
    nucl transcriptome
    """
    input:
        results = "blast_results/tblastn/c9g1i2_to_{taxon}_tblastn.txt",
        nucl_input_dummy = "txomes/{taxon}_nucl.nhr"
    output:
        best_hit = "blast_results/best_hit/tblastn/{taxon}_best_hit_tblastn.fasta"
    threads:
        1
    params:
        nucldb = "txomes/{taxon}_nucl"
    shell:
        """if [ -s {input.results} ]; \
        then blastdbcmd -db {params.nucldb} \
        -entry $(cat {input.results} | cut -f2 | head -1) > {output.best_hit} ;\
        else touch {output.best_hit}; \
        fi"""

rule tblastn_top_to_nr:
    """uses blastx to see if the top nucl hit from the transcriptome matches
    something in nr"""
    input:
        best_hit = "blast_results/best_hit/tblastn/{taxon}_best_hit_tblastn.fasta"
    output:
        results = "blast_results/blastx/{taxon}_tblastn_besttonr.txt"
    threads:
        maxthreads
    params:
        db = nrpath
    shell:
        """blastx -db {params.db} \
        -query {input.best_hit} \
        -query_gencode 1 \
        -outfmt 6 \
        -max_target_seqs 1 \
        -num_threads {threads} > {output.results}"""

rule get_top_blastx_hits_to_nr:
    """this gets the top hit, if it exists at all, from nr."""
    input:
        results = "blast_results/blastx/{taxon}_tblastn_besttonr.txt"
    output:
        best_hit = "blast_results/best_hit/{taxon}_tblastn_besttonr.fasta"
    threads:
        1
    params:
        db = nrpath
    shell:
        """if [ -s {input.results} ]; \
        then blastdbcmd -db {params.db} \
        -entry $(cat {input.results} | cut -f2 | head -1) > {output.best_hit} ;\
        else touch {output.best_hit}; \
        fi"""


# the next four rules deal with blasting the transcript of interest
#  against the transcriptomes of interest and extracting the sequences
rule blastp_c9g1i2:
    """this rule uses blastp to match the c9g1i2 transcript against
    the transcriptomes of interest. it outputs blast results to /blast_results
    """
    input:
        query = "fasta/luc_proteins/c9g1i2.pep",
        phr = "txomes/pep/{taxon}_pep.phr",
        pin = "txomes/pep/{taxon}_pep.pin",
        pog = "txomes/pep/{taxon}_pep.pog",
        psd = "txomes/pep/{taxon}_pep.psd",
        psi = "txomes/pep/{taxon}_pep.psi",
        psq = "txomes/pep/{taxon}_pep.psq"
    output:
        results = "blast_results/c9g1i2_to_{taxon}.txt"
    threads:
        maxthreads
    params:
        db = "txomes/pep/{taxon}_pep"
    shell:
        """blastp -db {params.db} \
        -query {input.query} \
        -outfmt 6 \
        -max_target_seqs 1 \
        -num_threads {threads} > {output.results}"""

rule get_top_hits:
    """ gets the peptide sequence of the top hit out of the target
    peptide transcriptome
    """
    input:
         results = "blast_results/c9g1i2_to_{taxon}.txt",
         pep_input_dummy = "txomes/pep/{taxon}_pep.phr"
    output:
        best_hit = "blast_results/best_hit/{taxon}_best_hit.pep"
    threads:
        1
    params:
        pepdb = "txomes/pep/{taxon}_pep"
    shell:
        """blastdbcmd -db {params.pepdb} \
        -entry $(cat {input.results} | cut -f2) > {output.best_hit}"""

rule blastp_top_to_nr:
    """blasts the peptide sequence of the top hit from the target transcriptome
    against nr"""
    input:
        best_hit = "blast_results/best_hit/{taxon}_best_hit.pep"
    output:
        results = "blast_results/{taxon}_best_to_nr.txt"
    threads:
        maxthreads
    params:
        db = nrpath
    shell:
        """blastp -db {params.db} \
        -query {input.best_hit} \
        -outfmt 6 \
        -max_target_seqs 1 \
        -num_threads {threads} > {output.results}"""

rule get_top_hits_to_nr:
    """this gets the top hit, if it exists at all, from nr."""
    input:
        results = "blast_results/{taxon}_best_to_nr.txt"
    output:
        best_hit = "blast_results/best_hit/{taxon}_best_to_nr.pep"
    threads:
        1
    params:
        db = nrpath
    shell:
        """if [ -s {input.results} ]; \
        then blastdbcmd -db {params.db} \
        -entry $(cat {input.results} | cut -f2 | head -1) > {output.best_hit} ;\
        else touch {output.best_hit}; \
        fi"""

rule orthology_report:
    """This prints out an orthology report of the species and its top
    hit to blast"""
    input:
        taxon_to_nr = sorted(expand("blast_results/{taxon}_best_to_nr.txt", taxon=config["species"])),
        taxon_to_nr_tblastn = sorted(expand("blast_results/blastx/{taxon}_tblastn_besttonr.txt", taxon=config["species"])),
        transcript_to_taxon = sorted(expand("blast_results/c9g1i2_to_{taxon}.txt", taxon=config["species"])),
        transcript_to_taxon_tblastn = sorted(expand("blast_results/tblastn/c9g1i2_to_{taxon}_tblastn.txt", taxon=config["species"]))

    output:
        orthoreport = "blast_results/orthology_report.txt"
    params:
        taxa = sorted(config["species"])
    run:
        outputf = open(output.orthoreport, "w")
        print("Species, blasttype, txid, Per_identity_w_c9g1i2, Aln_len_w_c9g1i2, E_value_with_c9g1i2, nr_id, Per_identity_w_nr, Aln_len_w_nr, E_value_with_nr", file = outputf)
        print(params.taxa)
        # do this once for blastp
        for blasttype in ["blastp", "tblastn"]:
            for i in range(len(params.taxa)):
                species = params.taxa[i]
                print("looking at {}".format(species))
                peridentx = ""
                alnlentx  = ""
                evaltx    = ""
                txid      = ""
                peridennr = ""
                alnlennr  = ""
                evalnr    = ""
                nrid      = ""
                if blasttype == "blastp":
                    inputfile = input.transcript_to_taxon[i]
                elif blasttype == "tblastn":
                    inputfile = input.transcript_to_taxon_tblastn[i]
                with open(inputfile) as f:
                    done = False
                    for line in f:
                        if line and not done:
                            line = line.split()
                            peridentx = line[2]
                            alnlentx  = line[3]
                            evaltx    = line[10]
                            txid      = line[1]
                            done = True
                if blasttype == "blastp":
                    inputfile = input.taxon_to_nr[i]
                elif blasttype == "tblastn":
                    inputfile = input.taxon_to_nr_tblastn[i]
                with open(inputfile) as f:
                    done = False
                    for line in f:
                        if line and not done:
                            line = line.split()
                            peridennr = line[2]
                            alnlennr  = line[3]
                            evalnr    = line[10]
                            nrid      = line[1]
                            done = True
                print("{}, {}, {}, {}, {}, {}, {}, {}, {}, {}".format(
                    species, blasttype,
                    txid, peridentx, alnlentx, evaltx,
                    nrid, peridennr, alnlennr, evalnr), file = outputf)
        outputf.close()

print(config["samples"])

rule hisat_build_to_transcript:
    input:
        reference = "fasta/isoform_analysis/c9g1i6.fasta"
    output:
        "alignments/c9g1i6_bt2index.1.ht2"
    threads:
        maxthreads
    params:
        indexbase = "alignments/c9g1i6_bt2index"
    shell:
        """
        hisat2-build -p {threads} {input.reference} {params.indexbase}
        """

rule map_rnaseq_unspliced_to_c9g1i6:
    input:
        reference = "fasta/isoform_analysis/c9g1i6.fasta",
        f_paired =   expand("reads/trimmed/{sample}_rna_f.trimmed.fastq.gz", sample = config["samples"]),
        r_paired =   expand("reads/trimmed/{sample}_rna_r.trimmed.fastq.gz", sample = config["samples"]),
    output:
        bamout = "alignments/{sample}_RNA_alignment_to_c9g1i6_bwa.sorted.bam",
        bamind = "alignments/{sample}_RNA_alignment_to_c9g1i6_bwa.sorted.bam.bai"
    threads:
        maxthreads
    shell:
        """bwa index {input.reference}; \

        bwa mem -t {threads} \
        {input.reference} \
        {input.f_paired} {input.r_paired} \
        | samtools view -b -F 4 -@ {threads} - \
        | samtools sort -@ {threads} - > {output.bamout}; \
        samtools index {output.bamout}"""

rule hisat_map_to_transcript:
    input:
        f_paired =   expand("reads/trimmed/{sample}_rna_f.trimmed.fastq.gz", sample = config["samples"]),
        r_paired =   expand("reads/trimmed/{sample}_rna_r.trimmed.fastq.gz", sample = config["samples"]),
        index_dummy_variable = expand("alignments/c9g1i6_bt2index.1.ht2", sample = config["samples"])

    output:
        bamout = "alignments/{sample}_RNA_alignment_to_c9g1i6.sorted.bam",
        bamind = "alignments/{sample}_RNA_alignment_to_c9g1i6.sorted.bam.bai"
    threads:
        maxthreads
    params:
        indexbase = "alignments/c9g1i6_bt2index"
    shell:
        """
        hisat2 -p {threads} -x {params.indexbase} \
        -1 {input.f_paired} -2 {input.r_paired} \
        | samtools view -h -b -F 4 -@ {threads} - \
        | samtools sort -@ {threads} - \
        > {output.bamout}; \

        samtools index {output.bamout}
        """

rule map_long_cdna_to_transcript:
    input:
        longcdna = [config["samples"][sample]["long_rna"] for sample in config["samples"]],
        reference = "fasta/isoform_analysis/c9g1i6.fasta"
    output:
        long_reads_bam = "alignments/{sample}_longcdna_to_c9g1i6.sorted.bam",
        long_reads_ind = "alignments/{sample}_longcdna_to_c9g1i6.sorted.bam.bai"

    shell:
        """minimap2 -t {threads} -ax splice {input.reference} {input.longcdna} \
        | samtools view -b -h -F 4 -@ {threads} - \
        | samtools sort -@ {threads} - \
        > {output.long_reads_bam};

        samtools index {output.long_reads_bam}
        """

rule hisat_map_pephits_to_transcript:
    input:
        fastq = "fasta/{sample}_pep_hits.fastq",
        reference = "fasta/isoform_analysis/c9g1i6.fasta"
    output:
        bamout = "alignments/{sample}_pepfrags_to_c9g1i6_bwa.sorted.bam",
        bamind = "alignments/{sample}_pepfrags_to_c9g1i6_bwa.sorted.bam.bai"
    threads:
        maxthreads
    params:
        alnname = "alignments/{sample}_c9g1i6"
    shell:
        """
        bwa index -p {params.alnname} {input.reference}; \
        bwa aln -t {threads} {params.alnname} {input.fastq} > {params.alnname}.sai; \
        bwa samse {params.alnname} {params.alnname}.sai {input.fastq} | samtools view -S -b -F 4 -@ {threads} - | samtools sort -@ {threads} - > {output.bamout}; \
        samtools index {output.bamout}
        """

rule hisat_map_pephits_to_final_transcripts:
    input:
        fastq = "fasta/{sample}_pep_hits.fastq",
        reference = "fasta/isoform_analysis/final_transcripts.fasta"
    output:
        bamout = "alignments/{sample}_pepfrags_to_final_transcripts_bwa.sorted.bam",
        bamind = "alignments/{sample}_pepfrags_to_final_transcripts_bwa.sorted.bam.bai"
    threads:
        maxthreads
    params:
        alnname = "alignments/{sample}_peps_to_final"
    shell:
        """
        bwa index -p {params.alnname} {input.reference}; \
        bwa aln -t {threads} {params.alnname} {input.fastq} > {params.alnname}.sai; \
        bwa samse {params.alnname} {params.alnname}.sai {input.fastq} | samtools view -S -b -F 4 -@ {threads} - | samtools sort -@ {threads} - > {output.bamout}; \
        samtools index {output.bamout}"""

rule trim_genomic_pairs:
    input:
        f1 = [config["samples"][sample]["genomic_f"] for sample in config["samples"]],
        f2 = [config["samples"][sample]["genomic_r"] for sample in config["samples"]],
        trim_jar = trimjar,
        adapter_path = adapterpath
    output:
        f_paired =   expand("reads/trimmed/{sample}_genomic_f.trimmed.fastq.gz", sample = config["samples"]),
        r_paired =   expand("reads/trimmed/{sample}_genomic_r.trimmed.fastq.gz", sample = config["samples"]),
        f_unpaired = temp(expand("reads/trimmed/{sample}_1.trim.unpaired.fastq.gz", sample = config["samples"])),
        r_unpaired = temp(expand("reads/trimmed/{sample}_2.trim.unpaired.fastq.gz", sample = config["samples"]))
    threads:
        maxthreads
    shell:
        """java -jar {input.trim_jar} PE \
        -threads {threads} \
        -phred33 {input.f1} {input.f2} \
        {output.f_paired} \
        {output.f_unpaired} \
        {output.r_paired} \
        {output.r_unpaired} \
        ILLUMINACLIP:{input.adapter_path}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36"""

rule trim_10X_forward:
    input:
        fastq = [config["samples"][sample]["10X_f"] for sample in config["samples"]],
        trim_jar = trimjar
    output:
        temp("reads/trimmed/{sample}_10X_f.temp.trim.fastq.gz")
    threads:
        maxthreads
    message:
        "Trimming the first 24 bases from {input.fastq}"
    shell:
        """java -jar {input.trim_jar} SE -threads {threads} \
        {input.fastq} {output} HEADCROP:24"""

rule trim_10X_reverse:
    input:
        fastq = [config["samples"][sample]["10X_r"] for sample in config["samples"]],
        trim_jar = trimjar
    output:
        temp("reads/trimmed/{sample}_10X_r.temp.trim.fastq.gz")
    threads:
        maxthreads
    message:
        "Trimming the first 3 bases from {input.fastq}"
    shell:
        """java -jar {input.trim_jar} SE -threads {threads} \
        {input.fastq} {output} HEADCROP:3"""

rule trim_10X_pairs:
    input:
        f1 = temp("reads/trimmed/{sample}_10X_f.temp.trim.fastq.gz"),
        f2 = temp("reads/trimmed/{sample}_10X_r.temp.trim.fastq.gz"),
        trim_jar = trimjar,
        adapter_path = adapterpath
    output:
        f_paired =   "reads/trimmed/{sample}_10X_f.trim.final.fastq.gz",
        r_paired =   "reads/trimmed/{sample}_10X_r.trim.final.fastq.gz",
        f_unpaired = temp("reads/trimmed/{sample}_10X_f.trim.unpaired.fastq.gz"),
        r_unpaired = temp("reads/trimmed/{sample}_10X_r.trim.unpaired.fastq.gz")
    threads:
        maxthreads
    message:
        "Trimming {input.f1} and {input.f2} as a pair."
    shell:
        """java -jar {input.trim_jar} PE \
        -threads {threads} \
        -phred33 {input.f1} {input.f2} \
        {output.f_paired} \
        {output.f_unpaired} \
        {output.r_paired} \
        {output.r_unpaired} \
        ILLUMINACLIP:{input.adapter_path}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36"""

rule trim_RNA_pairs:
    input:
        f1 = [config["samples"][sample]["rna_f"] for sample in config["samples"]],
        f2 = [config["samples"][sample]["rna_r"] for sample in config["samples"]],
        trim_jar = trimjar,
        adapter_path = adapterpath
    output:
        f_paired =   expand("reads/trimmed/{sample}_rna_f.trimmed.fastq.gz", sample = config["samples"]),
        r_paired =   expand("reads/trimmed/{sample}_rna_r.trimmed.fastq.gz", sample = config["samples"]),
        f_unpaired = temp(expand("reads/trimmed/{sample}_rna_1.trim.unpaired.fastq.gz", sample = config["samples"])),
        r_unpaired = temp(expand("reads/trimmed/{sample}_rna_2.trim.unpaired.fastq.gz", sample = config["samples"]))
    threads:
        maxthreads
    shell:
        """java -jar {input.trim_jar} PE \
        -threads {threads} \
        -phred33 {input.f1} {input.f2} \
        {output.f_paired} \
        {output.f_unpaired} \
        {output.r_paired} \
        {output.r_unpaired} \
        ILLUMINACLIP:{input.adapter_path}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 MINLEN:36"""
rule protpep_to_fasta:
    """this converts all of the peptide hits to a fasta file of sequences.
    This also makes a sam alignment file of hits to transcripts."""

    input:
        txome = [config["samples"][sample]["txome"] for sample in config["samples"]],
        protpep1 = "annotations/protein-peptides1.csv",
        protpep2 = "annotations/protein-peptides2.csv",
        protpep3 = "annotations/protein-peptides3.csv"
    output:
        fasta = "fasta/{sample}_pep_hits.fasta"
    run:
        hitindex = 0
        df1 = pd.read_csv(input.protpep1)
        df2 = pd.read_csv(input.protpep2)
        df3 = pd.read_csv(input.protpep3)
        dftemp = df1.append(df2, ignore_index=True)
        df = dftemp.append(df3, ignore_index = True)
        from Bio.Seq import Seq
        from Bio.SeqRecord import SeqRecord
        from Bio import SeqIO
        fastaout =  open(output.fasta, "w")
        indexer = 0
        seq_set = set()
        for record in SeqIO.parse(input.txome[0], "fasta"):
            print("         {}\r".format(indexer),end='')
            indexer += 1
            hits = df.loc[df["Protein Accession"] == record.id, ]
            if len(hits) > 0:
                for index, row in hits.iterrows():
                    start = row["Start"]
                    stop = row["End"]
                    hitstr = str(record.seq)[start:stop+1]
                    print(hitstr)
                    seq_set.add(hitstr)
        for uniqseq in seq_set:
            hitindex += 1
            newrecord = SeqRecord(Seq(uniqseq),
                                  id="hit{}".format(hitindex),
                                  description="")
            SeqIO.write(newrecord, fastaout, "fasta")
        fastaout.close()

rule fasta_to_fastq:
    input:
        expand("fasta/{sample}_pep_hits.fasta", sample = config["samples"])
    output:
        "fasta/{sample}_pep_hits.fastq"
    shell:
        """awk 'BEGIN {{RS = ">" ; FS = "\\n"}} NR > 1 {{print "@"$1"\\n"$2"\\n+"; \
        for(c=0;c<length($2);c++) printf "H"; printf "\\n"}}' \
        {input} > {output}
        """
