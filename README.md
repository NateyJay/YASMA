
# YASMA
### *<ins>Y</ins>et <ins>a</ins>nother <ins>sm</ins>all RNA <ins>a</ins>nnotator*
 
**warning**: this project is in active development and the readme may be out of date. 

### What is yasma?
This pipeline adds to a wide field of tools which have been used to assess ***small-RNA-sequencing*** data. Yasma is a **genome-based *de novo* approach** which is focused on **the whole sRNA population**, not just a few classes like miRNAs.

There are other approaches that follow a similar strategy, namely [ShortStack](https://github.com/MikeAxtell/ShortStack). Yasma try's to solve some persistent issues with these approaches, which appear to be exacerbated in challenging systems, such as sRNAs in Fungi.



### Problems with current approaches
* **Over-merging** of distinct, but closely oriented loci.
* **Creeping annotations** which don't don't represent the shape of a expressed region.
* **Under-merging** where numerous similar loci are annotated separately due to sequencing gaps.
* **Sensitivity to the depth** of a sRNA library relative to it's assembly size.


### General annotation strategy

Yasma relies on sRNA alignments to form genomic annotations. Alignments are performed using ShortStack3's alignment protocol (based on bowtie1).

Annotation based on these alignments follows a 3 step approach:
1. Building a ***sRNA coverage profile***.
2. Finding ***contiguous peaks*** of sufficient abundance.
3. Medium-distance merging of peaks which have ***similar sRNA profiles***.

This results in contiguous loci which are more homogenous in profile. It also tends to avoid over-annotation of background sequences.

These steps will be described in more detail below (*not yet implemented*).

# Yasma Modules

Yasma is broken into several modules, which form an analysis pipeline.

	yasma.py
	Usage: yasma.py [OPTIONS] COMMAND [ARGS]...

	YASMA (Yet Another small RNA Annotator)
	  
	Options:
	-h, --help  Show this message and exit.

	Commands:
	
	  Preliminary:
	    inputs      A tool to log inputs, which will be referenced by later tools.

	  Annotation:
	    peak        Annotator using peak identification and similarity merging.
	    
	  Calculation:
	    jbrowse     Tool to build coverage and config files for jbrowse2.
	    context     Compares annotations to identify cluster genomic context.
	    count       Gets counts for all readgroups, loci, strand, and sizes.
	    hairpin     Evaluates annotated loci for hairpin or miRNA structures.
	
	  Utilities:
	    adapter     Tool to check untrimmed-libraries for 3' adapter content.
	    readgroups  Convenience function to list readgroups in an alignment file.
***
#### Preliminary modules

***align*** is not implemented yet, but it will be a wrapper for shortstack's alignment functionality.

    yasma.py align -l lib_a.fa -l lib_b.fa -l lib_c.fa

***
***inputs*** is meant to log input file data, and do some basic prechecks. This function allows the user to input all data at this stage, and not have to call it again in subsequent functions.

    yasma.py inputs -o output_dir
    
***Result:*** 
This produces the `output directory` and `config.json`, logged with any values included as input. This file may be modified manually, following its format. Calls of following modules need only reference the `output_directory` to load these variables.
***
### Annotation module

***peak*** is the main annotation suite for yasma. This module will take sRNA alignments and form annotations, based on a user specified group of libraries (AKA readgroups).

	yasma.py peak -o output_dir -a alignment_file.bam -r wt_1 -r wt_2
	
***Result:***
This produces many files, representing the annotation. 
* `peak/loci.txt` and `peak/loci.gff3` - the result of the annotation, in tsv and gff formats. This is explained in detail below (*soon*).
* `peak/regions.gff3` - regions that were merged to produce loci.
* `peak/reads.txt` - reports the most expressed reads for a locus (up to a defined %).
* `peak/log.txt` - a log file of the stdout for the module.
* `peak/merges.txt` - a log of all of the merges that were performed.
* `peak/stats*` - run statistics for the chromosomes and overall.


### Calculation modules

Yasma has several modules which are used to extract key information  from the alignment result. These are broken into separate tools, to allow the user to analyze only those they wish to.
***
***context*** is used to identify the genomic context of sRNA loci, compared to an input gene_annotation. It used bedtools to find closest and intersecting features and returns a summary of them. 

	yasma.py context -o output_dir --gene_annotation_file ga_file.gff

***Result:***
`context/context.txt` is a tab delimited file summarizing the results of `bedtools closest` for loci against mRNAs, CDSs and exons. The result is summarized in the final column.

***
***count*** performs counting arithmetic on all libraries. 

	yasma.py count -o output_dir

***Result:***
This outputs 2 files: 
* `counts/counts.txt` - a count matrix of all loci by all libraries.
* `counts/deepcounts.txt` - a long format count table where abundances are further divided by sRNA-size and -strand -- useful for deeper analyses.
***
***hairpin***  tries to use folding with the sRNA profile to identify any sRNAs which might originate from hairpins. This includes miRNA testing using some basic published rules.

	yasma.py hairpin -o output_dir
	
***Result:***
* `hairpin/hairpins.txt` is a table summary of the findings of folding and miRNA identification attempts.
* `hairpin/folds/Cluster_*.eps` produces structural diagrams of the fold including expression profiles - mimicking the tool [strucVis](https://github.com/MikeAxtell/strucVis).
***
***jbrowse*** produces all key data for visualizing the annotation with jbrowse2. To properly view the output, this will need to be installed with it's base directory provided to the tool. When run it builds a complete configuration file which can be read by jbrowse2. If the folder is specified, it will add new entries to the`json`, and copy it with the data files them to the base folder.

	yasma.py jbrowse -o output_dir --jbrowse_directory base_directory/ -a alignment_file.bam -g genome_file.fa --gene_annotation_file ga_file.gff3

***Result:***
A full run of this module will produce `/jbrowse/config.json` - this config file, based on the specified jbrowse installation config.

Assembly and gene annotations from inputs are copied to a subfolder, named after the assembly.
* `/jbrowse/[assembly_name]/[assembly_name].fa`
* `/jbrowse/[assembly_name]/[assembly_name].fa.fai`
* `/jbrowse/[assembly_name]/[gene_annotation_file].gff3`

`loci` and `regions` annotations are added to sub_folder, based off of the annotation directory's name.
* /jbrowse/[assembly_name]/[output_direc

***
#### Utility modules
***adapter*** is a utility meant to help identify which are the most likely 3' adapter sequences in an untrimmed library. If a library is already trimmed, it will let the user know.

	yasma.py adapter -f library.fa

***
***readgroups*** is simply a wrapper for `samtools view -H` to show what libraries are found in 

	yasma.py readgroups -a alignment_file.bam
***
