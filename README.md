# SAMPhaser

```
Module:       SAMPhaser
Description:  Diploid chromosome phasing from SAMTools Pileup format.
Version:      0.8.0
Last Edit:    12/10/18
Copyright (C) 2016  Richard J. Edwards - See source code for GNU License Notice
```

## Function:

SAMPhaser is a tool designed to take an input of long read (e.g. PacBio) data mapped onto a genome assembly and
phase the data into haplotype blocks before "unzipping" the assembly into phased "haplotigs". Unphased regions are
also output as single "collapsed" haplotigs. This is designed for phasing PacBio assemblies of diploid organisms. 
By default, only SNPs are used for phasing, with indel polymorphisms being ignored. This is because indels are more 
likely to be errors. In particular, mononucleotide repeats could have indels that look like false well-supported polymorphisms.

### SAMPhaser overview

SAMPhaser first identifies variants from a pileup file generated using [SAMtools](https://github.com/samtools/samtools)
from a BAM file of mapped long reads. SNPs and indels are called for all positions where the minor allele is
supported by at least 10% of the reads (`mincut=X`), with an absolute minimum of two reads (`absmincut=X`).
The subset of biallelic SNPs with the minor variant supported by at least five reads (`absphasecut=X`) at a frequency
of at least 25% (`phasecut=X`) are used for phasing. Indels, and any SNPs not meeting these criteria, are used for
sequence correction, but not phasing.

Phasing is performed by iteratively assigning alleles and reads to haplotypes. Initially, each read is given an equal
probability of being in haplotype "A" or "B". The reference allele of the first SNP then defines haplotype A.
For each SNP, SAMPhaser iteratively calculates (1) the probability that each allele is in haplotype A given the
haplotype A probabilities for reads containing that allele, and then (2) the probability that each read is in
haplotype A given the haplotype A probabilities for that read's alleles at the last ten SNPs (`snpcalc=X`). This is
performed by modelling a SNP call error rate (`snperr=X` set at 5%) and then calculating the relative likelihood of
seeing the observed data if a read or allele is really in haplotype A versus haplotype B.

This progresses until all SNPs have been processed. If at any point, all reads with processed SNP positions reach
their ends before another SNP is reached, a new phasing block is started. Draft phase blocks are then resolved into
the final haplotype blocks by assigning reads and SNPs where the probability of assignment of a read to one haplotype
exceeds 95% (`trackprob=X`). Ambiguous reads and SNPs are ignored.

The final step is to "unzip" the reference sequence into "haplotigs". SAMPhaser unzips phase blocks with at least
five SNPs (`minsnp=X`). Regions that are not unzipped are output as "collapsed" haplotigs. First, phased reads are
assigned to the appropriate haplotig. Regions of 100+ base pairs without coverage (`splitzero=X`) are removed as
putative structural variants, and the haplotig split at this point. Haplotigs with an average depth of coverage
below 5X (`minhapx=X`) are removed. Note that this can result in "orphan" haplotigs, where the minor haplotig did not
have sufficient coverage for retention. Haplotigs ending within 10 bp (`endmargin=X`) of the end of the reference
sequence are extended. Next, collapsed blocks are established by identifying reads that (a) have not been assigned to
a haplotype, and (b) are not wholly overlapping a phased block.

Finally, unzipped blocks have their sequences corrected. This is performed by starting with the reference sequence 
and then identifying the dominant haplotype allele (or consensus for collapsed blocks) at all variant positions 
(not just those used for phasing) providing the variant has at least 10% (min. three) reads supporting it 
(`unzipcut=X` `absunzipcut=X`). The final haplotig sequence is the original reference sequence with any assigned 
non-reference alleles substituted in at the appropriate positions. Single base deletions are cut out of the sequence 
and so it may end up shorter than the original contig. Insertions and longer deletions are not currently handled and 
are ignored; for this reason, it is important to re-map reads and correct the final haplotig sequences.

## Running SAMPhaser

To install, simply download or clone either this repository or the main [SLiMSuite](https://github.com/slimsuite/SLiMSuite) repository. SAMPhaser is written in Python 2.x and can be run directly from the commandline:

    python $CODEPATH/samphaser.py [OPTIONS]

If running as part of SLiMSuite, `$CODEPATH` will be the SLiMSuite `tools/` directory. If running from the standalone SAMPhaser git repo, `$CODEPATH` will be the path the to `code/` directory. 

The basic SAMPhaser run command needs a genome sequence (`seqin=FASFILE`) and pileup file (`pileup=FILE`):

    python $CODEPATH/samphaser.py -seqin <genome.fasta> -pileup <genome.pileup>

## Documentation

Documentation is available in the [SAMPhaser.md](./SAMPhaser.md) file included in this repository. A list of commandline options can also be generated by running with the `-help` option.

## Citing SAMPhaser

SAMPhaser is not yet published. If you want to use SAMPhaser in a publication in the meantime, please cite the main [SLiMSuite release](https://github.com/slimsuite/SLiMSuite/releases) Zenodo [DOI](https://zenodo.org/record/1302990).

