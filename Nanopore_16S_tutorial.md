# Nanopore 16S tutorial

Current available tools for nanopore based 16S sequence are currently not up to short reads standard. For instance you would be hard pressed to obtain asv level resolution. It is still however quite easy to derive a genus level characterization of you samples.

In this tutorial we will show you how to use a reference based approach to generate taxonomic profiles of your nanopore 16S samples

The workflow is quite typical and involve

1. [Read mapping](#readmapping)

2. [alignment filtering](#al)

3. [Phyloseq in R ](#phyloseq)
 
## Getting started (env)


We use a [conda](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf) env to install all dependencies, you don't need to install anything and all dependencies are available but only inside that environment.   

Try to check the command line help of minimap2

	    minimap2 -h
<details><summary>not working?</summary>
<p>
Conda environment are created as independant environment to everything else, you need to "activate" an environment to be able to access the sets of tools installed inside.

	    conda env list
	    conda activate workshop
	    minimap2 -h

</p>
</details>

## Reads mapping
Let's create a Projects directory and work inside:
```bash
	mkdir $HOME2/Mystery_16S
	cd $HOME2/Mystery_16S
```
We use minimap2 to map reads to the silva database. As with most bioinformatic tools, information on how to run this tools can be obtain either by google or by typing `minimap2 -h`
Please spend at least 5 minute trying to run minimap2 with as read file:
```bash 
	$DATA/mystery_16S/BC01.fastq.fq
```
and using silva database found at:
```bash 
	/home/ubuntu/software/silva_db/No_U_SILVA_138.fa
```
How to choose the preset? 
<details><summary>the correct command is:</summary>
<p>

```bash
	cd $HOME2/Mystery_16S
	mkdir Map
	minimap2 -a -x map-ont  /home/ubuntu/software/silva_db/No_U_SILVA_138.fa $DATA/mystery_16S/BC01.fastq.fq >Map/BC01.sam
```
</p>
</details>

You can look at the sam:
```
	tail Map/BC01.sam
```

It is quite a complex [format](https://en.wikipedia.org/wiki/SAM_(file_format))

The sam file is a bit bulky so we never store alignments in this format instead we would convert it into bam. Can you convert this file using 
'samtools view':


<details><summary> Convert sam to bam command</summary>
<p>

```bash
    cd Map
    samtools view -h -b -S BC01.sam > BC01.bam
```
</p>
</details>

Using samtool we can filter only those reads which are mapped to the silva database.
```bash
    samtools view -b -F 4 BC01.bam > BC01.mapped.bam
```

For downstream analysis we needs the bam file to be sorted and indexed
```bash
	samtools sort BC01.mapped.bam -o BC01.mapped.sorted.bam 
	samtools index BC01.mapped.sorted.bam 
```
What is indexing/sorting? 

To run all samples we would place these steps in a shell script:

```bash
	cd $HOME2/Mystery_16S
	rm $HOME2/Mystery_16S/Map/*

	for file in $DATA/mystery_16S/*.fq
	do 
		name=${file%.fastq.fq}
		name=${file#*mystery_16S*}		
		minimap2 -a -x map-ont  /home/ubuntu/software/silva_db/No_U_SILVA_138.fa $file | samtools view -b -F 4 - | samtools sort - > Map/$name.mapped.sorted.bam
		samtools index Map/$name.mapped.sorted.bam
	done
```
## Alignment filtering
This next step is relatively quite simple as this is a in house script which does multiple step in a unique command line.
It uses some python library to interface with the .bam file, select best hits and filter them. It also gather taxonomic information from silva database and format it in a systematic way. 
Finally it output 2 tables which can be used for phyloseq type analysis, namely a taxonomy file which lists all taxa observed and a count file with the corresponding counts. 
To run it, simply ..... try to make sense of the help message

<details><summary> Do try a bit before uncovering</summary>
<p>

```bash
	cd $HOME2/Mystery_16S
	mkdir phyloseq
	/home/ubuntu/software/Respharm_ont_workshop/map_taxa.py Map /home/ubuntu/software/silva_db/No_U_SILVA_138.fa phyloseq 
```

</p>
</details>

What is a script? This is a python script how do you look inside? What are the library used?


## Phyloseq in R
For the next step we are going to use R and in particular the library [phyloseq](https://joey711.github.io/phyloseq/index.html). Multiple more complete tutoriels exists on this website and [this one](https://joey711.github.io/phyloseq/import-data.html) is a good one to start with. 

First let's load the library needed
```
	library(phyloseq)
	library(data.table)
	library(ggplot2)
```

Lets load in R the 2 tables we just created. Please replace "PATH" by the correct value. 

    taxa = read.table("PATH",header=TRUE,sep="\t",na.strings="*",stringsAsFactors=F,row.names=1,check.names=FALSE)
    taxa[taxa==""]="NA"
    cnts = read.table("PATH",header=TRUE,sep="\t",na.strings="*",stringsAsFactors=F,row.names=1,check.names=FALSE)

There is quite a number of option there, what do they mean? 

Use in turn the following command on both tables:

    dim(taxa)
    str(taxa)
    head(taxa)

How many different taxa did we detect?

Some of theses are quite low frequency and are quite likely noise/errors. We are removing anything whith less than 20 counts over all samples.

    selection = rowSums(cnts)>10
    taxa = taxa[selec,]
	cnts = cnts[selec,]

What is the new number of taxa? 

We will now create a phyloseq object. 

    taxa_ps = tax_table(taxa)
    colnames(taxa_ps)=colnames(taxa)
    
    cnt_tbl = otu_table(t(cnts), taxa_are_rows=FALSE)
    ps = phyloseq(cnt_tbl,taxa_ps)

If you type ps what do you see?
The most difficult part is now done, we can now use all of phyloseq functions and realize multiple plots.

    pdf("taxa_bar_plots.pdf")
    ps_king = tax_glom(ps, "kingdom")
    plot_bar(ps_king, fill ="kingdom")
    ps_P = tax_glom(ps, "phylum")
    plot_bar(ps_P, fill="phylum")
    ps_C = tax_glom(ps, "class")
    plot_bar(ps_C, fill="class")
    ps_O = tax_glom(ps, "order")
    plot_bar(ps_O, fill="order")
    ps_F = tax_glom(ps, "family")
    plot_bar(ps_F, fill="family")
    ps_G = tax_glom(ps, "genus")
    plot_bar(ps_G, fill="genus")
    dev.off()

To see this pdf you will need to use either Winscp/cyberduck on windows or just scp on mac/ubuntu.
