# annotator
 

A general look at better ways to define loci in sRNA seq data.


### Problems with current approaches
* ShortStack does a good job, with a easily understood algorithm and set of rules. However, it misses lots of miRNAs because of concatenation of loci.
* This problem is exacerbated in fungi - high background sRNA levels make it hard to distinguish off-sized, degradation-related, and other loci from dicer or RNAi-related loci.
* Moreover, ShortStack defines the boundaries of loci based on a merge window that require only a single read to extend. This means in libraries which are overly-deep compared to the genome size, loci can be come trivially large (annotating the whole genome). This is what we see in several fungi.
* **The problems we see in fungi should be the emphasis of this software**


### Ideas for implementation
*Some implemented, others not...*

#### Loci and sub-loci.
* one approach could be to try to only provide more detail to those defined by shortstack. Here, we would define loci which are within shortstack loci, hopefully finding shorter loci which are miRNAs.
* This might be a good way to find miRNAs, but I think it would still be very problematic for fungi. Looking at the Tratr/Bocin annotations from Consuelo's work - virtually the entire genome is loci.


#### Expression thresholding
* Perhaps we could have a high minimum locus coverage threshold to find a locus, with less strict thresholds for it's edges. 
* I don't know if this would directly address the problems shortstack has with miRNA discovery.
* Binning (not this actually) or a density-driven algorithm could be a basis for this.
* A quantile-based RPM threshold?


#### miRNA-specific algorithms
* Highly-multiplexed folding would allow us to try to find resilient hairpins within loci. *Are these representative of real hairpins?* It looks like from the fungi hpRNA paper that we can recover miRNAs pretty well in control species using this approach. In any case, I think more folding is an answer to this problem.
* Maybe we could fold stretches that are highly stranded with multiple sizes. Most frequent hits could be condensed to find the most resilient hairpin.



## Current approach

The current version of the code is generally similar to ShortStack, but implements some of the changes above and some more suitable defaults to account for the changes. The pipeline currently follows 3 steps:

#### **1)** Alignment


#### **2)** Precheck

#### **3)** Annotation









