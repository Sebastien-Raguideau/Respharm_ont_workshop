# nanopore_genomic_tutorial 

This tutorial is *not* 100% step-by-step, since you will have to navigate to the right folders, make directories, maybe rename a file or two, etc. When trying to run a tool, first check the link given for usage instructions or use the --help flag (e.g. `flye --help`) to get instructions. I have hidden some commands - the idea is that you can try and figure out the settings yourself by checking the documentation. The solutions can be seen by clicking on the arrow

<details>
<summary>Like here</summary>
       
    scp --help
    
</details>

The nanopore outputs are in ~/nanopore_data. fastq_pass is the directory containing the basecalled reads meeting the selection criteria. The report file contains some good info on the sequencing run. In ~/outputs you will find outputs from all tools. All samples have been processed with filtlong, assembled with flye and polished with medaka, but all other tools have only been applied to one or two samples.

If you're unsure about the file formats fasta, fastq, sam & bam, check out: https://bioinformatics.uconn.edu/resources-and-events/tutorials-2/file-formats-tutorial/#

Everything needed for the workshop is installed using conda. Conda makes it very easy to install and run programs. For running a program, its environment has to be activated, like so:

    conda activate flye

And you'll see the (flye) environment being active in the line where you type.

When you're done with the program, deactivate with

    conda deactivate

To see a list of all environments on the server, type

    conda env list

For this workshop, run all commands with a maximum of 6 threads.

## Read QC
Long reads are *already* basecalled and adapters trimmed (done with guppy basecaller). The raw basecalled reads are located in `~/nanopore_data/fastq_pass/barcode[18/19/20/21/22]/`. The read files were put into one file (e.g `barcode18.fastq`) using the cat command (`cat ./*.fastq > barcode_xy_all_reads.fastq`). For this example, all reads _are already in one file_

Let's check out the quality of the long reads before assembly.

To check for possible contamination in our sample, you can use kraken2 to classify all the reads & provide a report, then check report - how many % are E. coli? Any unexpected things? https://ccb.jhu.edu/software/kraken2/

<details>
<summary>Kraken command</summary>
    
    kraken2 --db ~/software/kraken_db/ --threads 6 --output barcode_18.output --report barcode_18.report --use-names ~/nanopore_data/fastq_pass/barcode18.fastq

</details>

    
You can also visualise the reports using pavian in Rstudio on your laptopn if you have it. Not super necessary for this workshop.

### Seqkit

Use seqkit stats to check stats on N50, number of reads, max length, quality scores, etc. Compare the different fastq files - notice any differences? Which barcode has longer reads, which one has the most sequences? What is the estimated depth of coverage for each sample when taking into account E. coli genome size? https://bioinf.shenwei.me/seqkit/usage/#stats

    seqkit stats -a barcode18.fastq

### NanoPlot (optional)
Use NanoPlot to make nice graphs with read length and quality, see usage here: https://github.com/wdecoster/NanoPlot
You need to download the graphs using scp to look at them. Do the graphs agree with the info you got from seqkit?

<details>
<summary>NanoPlot command</summary>
    
    NanoPlot --fastq ~/nanopore_data/fastq_pass/barcode19.fastq --loglength -o barcode19 --threads 6

</details>

### Filtlong
Some of the barcodes have extremely high coverage! Assembling a >200x coverage genome takes about 45 minutes. Reducing it to 100 makes it much quicker, only ca. 20 minutes. You can use filtlong to remove the short and low-quality reads, thereby reducing assembly time. Check the documentation: https://github.com/rrwick/Filtlong
One ting to take into account - how many bases should we keep to have 100x coverage?

<details>
<summary>FiltLong command</summary>
    
    filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 barcode22.fastq | gzip > barcode22_filtered.fastq.gz

</details>


## Assemblies
We're going to do  long-read assembly followed by polishing with medaka

## Long Read Assembly with Flye
Assemble long reads with flye. This usually produces the most contiguous assembly, but it still needs polishing both with long-read specific tools and with the short reads of the same sample. This should take 20 minutes or so. Type flye -h to get instructions, or check https://github.com/fenderglass/Flye

<details>
<summary>Flye command</summary>
       
    flye --threads 6 --genome-size 4m --nano-raw barcode18.fastq --out-dir barcode18
    
</details>

You can inspect the assembly outcome by checking assembly_info.txt. Is there a circular chromosome? Any plasmids? Any assembly artifacts?

    less flye_assembly/assembly_info.txt


## Long Read Assembly Polishing with Medaka

Do not polish with racon - _new versions of medaka now work with the output of flye directly._


### Medaka
Medaka is a machine-learning model based polisher, so it is _super important_ to correctly specific the flow cell and the basecaller used, because that will influence the way that medaka corrects the assembly. You can find the flow cell, the guppy basecaller version and the basecalling model used in the sequencing report file in the nanopore data folder. The available models can be found by typing medaka_consensus -h
The main output of medaka consensus.fasta in the medaka folder.
https://github.com/nanoporetech/medaka#Usage

<details>
<summary>Medaka command</summary>
    
    medaka_consensus -m r941_min_hac_g507 -t 12 -i ~/outputs/filtlong/barcode18.fastq_filtered.fastq.gz -d ~/outputs/flye/barcode18/assembly.fasta -o barcode18_medaka
    
</details>

## Assessing the assembly
When we have an assembly, we should check how complete it is. That can be done in different ways.

### Bandage
Bandage visualises the assembly graph. It it similar info as you will get from assembly_info.txt, but including nice visualisations. Download the assembly_graph.gfa file from the flye output onto your computer and then load into Bandage: https://github.com/rrwick/Bandage/

### Seqkit
Seqkit can also give you some statistics on your assembly:

       seqkit stats -a consensus.fasta

### BUSCO  
Busco uses single copy core genes (SCCGs) to assess completeness and assembly quality. We can use BUSCO to assess the different assemblies we created, as well as to show the improvement of the polishing steps in the course of the nanopore assembly. https://busco.ezlab.org/

First, select the appropriate dataset listed by BUSCO that fits with the phylogeny of our *Knoellia* genome. We should use the most specific dataset we can, so we can list all available datasets and choose the most specific one:

    busco --list-datasets
    
Another option is to use the --auto-lineage flag, then BUSCO will try selecting the best lineage itself. Alternatives to BUSCO, especially for MAGs, include checkm (more computing power required).

The result will show how many of the expected SCCGs are present & complete (C), if they are fragmented (F), or if they are missing completely (M). This indicates whether a genome is complete or a part is missing; it also gives a measure of sequence quality, since lots of fragmented SCCGs mean high rates of indels.

Comparing different assemblies and taking into account what you know about the read quality, what are your conclusions?

<details>
<summary>BUSCO command</summary>    
    
    busco -m genome -i ~/outputs/flye/barcode18/assembly.fasta -o barcode18_flye_only ---auto-lineage  -c 6
    busco -m genome -i ~/outputs/medaka/barcode18_medaka/consensus.fasta -o barcode18_medaka --auto-lineage -c 6
</details>
    
### IDEEL
Ideel is a way of assessing the prevalence of frameshift mutations resulting from insertion/deletion errors in your assembly. It takes a while to run! So maybe run it later or on the side. https://github.com/mw55309/ideel

To run it:
Put your genomes in a directory with the name "genomes" in the ~/programs/ideel directory, makes sure they have ".fa" file extension. Then navigate to the ideel directory, activate conda environment and execute

    snakemake


