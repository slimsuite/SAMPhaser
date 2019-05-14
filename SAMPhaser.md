# SAMPhaser Documentation 

This is the draft documentation for:

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

## Commandline:

    ### ~ Input/Output options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    seqin=FASFILE   : Input genome to phase variants in []
    pileup=FILE     : Pileup file of reads against input genome. []
    basefile=FILE   : Root of output file names (same as Pileup input file by default) []
    mincut=X        : Minimum read count for minor allele (proportion of QN if <1) for pileup parsing [0.1]
    absmincut=X     : Absolute minimum read count for minor allele (used if mincut<1) [2]
    indels=T/F      : Whether to include indels in "SNP" parsing [True]
    snptableout=T/F : Output filtered alleles to SNP Table [False]
    skiploci=LIST   : Optional list of loci (full names or accnum) to skip phasing []
    phaseloci=LIST  : Optional list of loci (full names or accnum) to phase (will skip rest) []
    ### ~ Phasing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    phasecut=X      : Minimum read count for minor allele for phasing (proportion of QN if <1) [0.25]
    absphasecut=X   : Absolute minimum read count for phasecut (used if phasecut<1) [5]
    phaseindels=T/F : Whether to include indels in "SNP" phasing [False]
    snperr=X        : Probability of an incorrect (biallelic) SNP call for individual read nucleotides [0.05]
    snpcalc=X       : Max number of SNPs to use for read probability calculations (fewer = quicker) [10]
    trackprob=X     : Min probability for assigning a read/SNP to Track A/B [0.95]
    ### ~ Unzipping options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###
    minsnp=X        : Min number of SNPs per phased haplotype block [5]
    endmargin=X     : Extend block ends within X nucleotides of sequence ends [10]
    unzipcut=X      : Minimum read count for allele for unzipping haplotigs (proportion of QN if <1) [0.1]
    absunzipcut=X   : Absolute minimum read count for unzipcut (used if unzipcut<1) [3]
    minhapx=X       : Minimum mean coverage for haplotig [5]
    halfhap=T/F     : Whether to allow "half haplotigs" where one halpotig in a pair is removed by minhapx [True]
    splitzero=X     : Whether to split haplotigs at zero-coverage regions of X+ bp (-1 = no split) [100]
    rgraphics=T/F   : Whether to generate PNG graphics using R. (Needs R installed and setup) [True]
    poordepth=T/F   : Whether to include reads with poor track probability in haplotig depth plots (random track) [False]
    ### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ###

## History
    # 0.0.0 - Initial Compilation.
    # 0.1.0 - Updated SAMPhaser to be more memory efficient.
    # 0.2.0 - Added reading of sequence and generation of SNP-altered haplotype blocks.
    # 0.2.1 - Fixed bug in which zero-phasing sequences were being excluded from blocks output.
    # 0.3.0 - Made a new unzip process.
    # 0.4.0 - Added RGraphics for unzip.
    # 0.4.1 - Fixed MeanX bug in devUnzip.
    # 0.4.2 - Made phaseindels=F by default: mononucleotide indel errors will probably add phasing noise. Fixed basefile R bug.
    # 0.4.3 - Fixed bug introduced by adding depthplot code. Fixed phaseindels bug. (Wasn't working!)
    # 0.4.4 - Modified mincut=X to adjust for samtools V1.12.0.
    # 0.4.5 - Updated for modified RJE_SAMTools output.
    # 0.4.6 - splitzero=X : Whether to split haplotigs at zero-coverage regions of X+ bp (-1 = no split) [100]
    # 0.5.0 - snptable=T/F    : Output filtered alleles to SNP Table [False]
    # 0.6.0 - Converted haplotig naming to be consistent for PAGSAT generation. Updated for rje_samtools v1.21.1.
    # 0.7.0 - Added skiploci=LIST and phaseloci=LIST  : Optional list of loci to skip phasing []
    # 0.8.0 - poordepth=T/F   : Whether to include reads with poor track probability in haplotig depth plots (random track) [False]
    
## How SAMPhaser works

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

**NOTE:** This documentation is in development. Please report any anomalies.

## Running SAMPhaser

The basic SAMPhaser run command needs a genome sequence (`seqin=FASFILE`) and pileup file (`pileup=FILE`):

    python ~/slimsuite/tools/samphaser.py -seqin <genome.fasta> -pileup <genome.pileup>

SAMPhaser then runs 4 main phases:

1. [Pileup parsing](#parsing).
2. [Phasing](#phasing).
3. [Unzipping](#unzipping).
4. [Report generation](#report).

**NOTE:** If the `*.blocks.tdt` output file is found, SAMPhaser will load this and skip straight to the unzipping step. Otherwise, if the `rje_samtools` output (below) is already present, it will be read in and reused (skipping the parsing step) unless `force=T`.


<a name="input"></a>
### SAMPhaser input

**Genome fasta file.** This should be in standard fasta format. Most naming formats should work, but use [SLiMSuite format](http://slimsuite.blogspot.com/2015/10/file-format-fasta-seqfile-fasfile.html) for best results.

**Pileup file.** This should have reads mapped again the genome file given. Example commands to generate this file are given in the **Pileup parsing** section.

If the `*.blocks.tdt` output file is found, SAMPhaser will load this and skip straight to the unzipping step. Otherwise, if the `rje_samtools` output (below) is already present, it will be read in and reused unless `force=T`.

<a name="output"></a>
### SAMPhaser output

Pileup parsing by `rje_samtools` generates the following outputs:

* `*.rid.tdt` = A table of individual read start/stop positions. (**NOTE:** SAMPhaser expects this to have been generated from the sorted pileup file. If you have previously generated this from an unsorted SAM file, you may need to delete it so that it is remad.

            RID	Locus	Start	End
            176	chrI_YEAST__BK006935	61	351
            ...

* `*.Q30.10.snponly.bFiT.tdt` = A table of SNP frequencies. This suffix is generated by default settings. The `QX.N.` part is from the `rje_samtools` settings `qcut=X` and `minqn=X`. `bXiX` marks whether `biallelic=T/F` and `indels=T/F`.

            Locus	Pos	Ref	N	QN	Seq	Dep	RID
            chrI_YEAST__BK006935	26	A	148	148	A:107|C:2	0.74	A:2,3,6,7,8,9,10,12,14,15,17,18,19,20,21,22,23,25,28,29,30,32,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,68,70,71,73,75,76,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,98,99,100,102,103,104,106,108,110,112,116,117,119,120,121,123,124,125,128,129,131,132,135,136,137,138,140,142,146,147,148|C:78,127
            ...

SAMPhaser then generates the following additional output files. By default, file prefixes match that of `pileup=FILE`. This is over-ridden with `basefile=X`, which will also set the logfile name unless `log=X` is also used.

* `*.snp.tdt` = A table of the phased SNPs and their read/block assignments.

        Locus	Pos	Ref	N	QN	Seq	Dep	V1	R1	A1	V2	R2	A2	Block

* `*.haprid.tdt` = A table of phased reads and their SNPs.

       RID	Locus	Start	End	Block	Track	pTrack	SNP

* `*.hapsnp.tdt` = An updated table of the phased SNPs and their phasing.

       Locus	Block	Pos	Ref	A	B	nA	nB	pA

* `*.blocks.tdt` = A table of the parsed haplotype blocks and the SNP calls for each one.

       Locus	Block	Track	Start	End	HapStart	HapEnd	SNP

* `*.depthplot.tdt` = Read depths across the assembly, for plotting.

       Locus	Pos	X

* `*.coverage.tdt` = Coverage statistics for each locus (scaffold/contig).

       Locus	Length	MeanX	MinX	MaxX	MedianX	Coverage

* `*.haplotigs.tdt` = A table of each haplotig, with position, SNP and read counts.

        Locus	Block	Track	Start	End	HapStart	HapEnd	Haplotig	HapAcc	nSNP	nRID	MeanX

* `*.haplotigs.fas` = Fasta files of the phased and unzipped haplotigs.
* `*.haplotigs.rid.tdt` = Table of reads assigned to each haplotig.

        RID	Locus	Start	End	pTrack	Haplotig	HapAcc

* `*.haplotigs.depthplot.tdt` = Read depths across the haplotigs, for plotting.

        Locus	Pos	X

* `*.haplotigs.coverage.tdt` = Coverage statistics for each haplotig.

      Locus	Length	MeanX	MinX	MaxX	MedianX	Coverage

* `*.samphaser.html` = Final SAMPhaser report. 
* `*.samphaser.log` = Log of SAMPhaser processing.

_**NOTE:** SAMPhaser is still under development. Example data and improved documentation will follow. Please contact the author in the meantime with questions._

<a name="parsing"></a>
### Pileup Parsing
SAMPhaser uses the `rje_samtools.SAMTools` class (`rje_samtools.py`) to parse a `*.pileup` file created from a SAM file by SAMTools (example with `seqin=genome.fasta pileup=genome.pileup`):

```
samtools faidx genome.fasta                          # Index genome
samtools view -bo genome.bam -S genome.sam"          # Convert SAM to BAM
samtools sort -@ 16 -f genome.bam genome.sort.bam    # Sort BAM file (16 threads)
samtools index -b genome.sort.bam                    # Index sorted BAM file
samtools mpileup -BQ0 -d10000000 -f genome.fasta genome.sort.bam > genome.pileup  # Make Pileup file.
```

`rje_samtools` is given a new default setting of `mincut=0.1`. This means that the default parsed file will contain all variable positions where the minor allele frequency is at least 10% (`mincut=X`) and has at least 2 reads (`absmincut=X`). SAMPhaser needs the `*.rid.tdt` file generated by `rje_samtools` and thus sets: `rid=T snponly=T`. 

Although SAMPhaser by default only considers biallelic SNPs, `biallelic=F` is set by default as non-biallelic positions are used during the unzipping process to modify the haplotig sequences. Likewise, `indels=T` is used by default, which affects the unzipping process but not phasing. The `phaseindels=T/F` controls whether to include indels in "SNP" phasing, [default: False]. **NOTE:** `indels=F` will also set `phaseindels=F` by removing indels from the parsed pileup file.

Following `rje_samtools` parsing, SAMPhaser performs additional filtering of the parsed SNPs, reducing them to postions with _exactly 2 alleles_ with a minimum frequency of `phasecut=X` (default=0.25) and minimum read count of `absphasecut=X` (default=5). If `phaseindels=F` (default) then positions are restricted to those with two non-gap characters.


<a name="phasing"></a>
### SNP Phasing

Following SNP parsing, SAMPhaser performs the phasing step. This is done by iteratively assigning alleles and reads to haplotypes, which are labelled `A` and `B`.

The SNP table has the following fields added:

* `V1` = Allele 1
* `R1` = List of RID for allele 1
* `A1` = Probability of Track A for Allele 1
* `V2` = Allele 2
* `R2` = List of RID for allele 2
* `A2` = Probability of Track A for Allele 2
* `Block` = Assigned haplotype block

The RID table has the following fields added:

* `A` = Prob of trackA
* `B` = Prob of trackB
* `SNP` = List of SNPs (pos,allele)
* `Block` = Haplotype Block

Initially, each read is given an equal 50% probability of being in track A and B and alleles of first SNP are assigned to tracks. The reference allele of the first SNP defines track A.

SAMPhaser works through each locus in turn and then through each SNP position in order. For each position, SAMPhaser iteratively calculates the probability that each allele is in Track A given the Track A probabilities for reads containing that allele, and then the probability that each read is in Track A given the Track A probabilities for the alleles that read contains. (See `snpA` and `readA` calculations, below.

Each time the read probabilities are calculated, the last `snpcalc=X` SNPs have their allele assignment probabilities recalculated. This also sets the number of alleles used to calculate read assignment probabilities. By default, the last 10 SNPs are used for this purpose.

_Increasing `snpcalc=X` is expected to make calculations more accurate at the cost of speed. Future work may compare different `snpcalc=X` settings more directly._

This progresses until all SNPs have been processed. If at any point, all reads with 1+ processed SNP positions reach their ends before another SNP is reached, a new phasing block is started.

This ends with the log statement for each LOCUS:

    #PHASE - Phased X LOCUS SNPs: N phased total haplotype blocks

The next step is to resolve these draft blocks into the final haplotype blocks where there is sufficient support. Where the probability of assignment of a read to TrackA or TrackB meets the `trackprob=X` cutoff (95% by default), that read is assigned to that block. SNPs are likewise assigned to tracks if the probability satisfies `trackprob=X`. The extent of the block is then calculated. This corresponds to two distinct length concepts:

1. The phase block itself is defined by the extreme positions of the SNPs in that block (`HapStart` and `HapEnd`). This will be identical for each Track.

2. The track block is defined by the extreme starts and ends of the mapped reads assigned to that track (`Start` and `End`). Each track will extend beyond each end of the phase block. Adjacent blocks will generally overlap in terms of read coverage and thus have overlapping ends, even though a lack of SNPs between them allows the phase blocks to merge.

Phasing information is output to `*.haprid.tdt` and `*.hapsnp.tdt` (see above).


**snpA calculation.** This sets the probability of each allele in the SNP position under consideration being trackA, using:

1. The list of RIDs currently being phased.
2. The probability of an incorrect (biallelic) SNP call for individual read nucleotides, set by `snperr=X` (default=0.05).

If there are no RIDs currently being phased, a new phase block will be started. In this case, the Reference allele (_e.g._ the one matching the `seqin=FILE` sequence) is defined as Track A. If neither allele matches the reference, Allele 1 is set as Track A.

Assuming there are reads currently being phased, the raw TrackA:TrackB likelihood ratio is calculated for each allele. For each allele:

* `pA` = exp( sum-over-reads-with-allele[ math.log(rp*(1-erate) + (1-rp)*erate) ] )
* `pB` = exp( sum-over-reads-with-allele[ math.log((1-rp)*(1-erate) + rp*erate) ] )

where `rp` is the probability that a given read is Track A and erate is the probability of an incorrect (biallelic) SNP call set by `snperr=X`. For each calculation, only reads with the alleles under consideration at that position are included, _i.e._ sequencing errors will be ignored from the calculation (unless they generate the _wrong_ allele).

This works on the principle that each read has an `erate` probability of having the "wrong" allele. Therefore the TrackA probability for a given read is the combination of the TrackA read probability if the allele is correctly assigned, `rp*(1-erate)`, plus the probability that the read is TrackB and the allele is incorrect, `(1-rp)*erate`. The overall likelihoods are the product of these values over all reads with that allele.

The likelihood ratios of `pA` and `pB` are then normalised to sum to 1.0. The TrackA (`pa`) probabilities for the two alleles are then normalised to sum to 1.0 and the normalised value assigned as the `snpA` probability for that allele.


**readA calculation.** This sets the probability that a read is in Track A, given:

1. The list of alleles in the read. This is restricted to the last `snpcalc=X` SNP positions that have been calculated.
2. The probability of an incorrect (biallelic) SNP call for individual read nucleotides, set by `snperr=X` (default=0.05).

For each assessed SNP in the read (up to the last `snpcalc=X`), the probability of that allele being from each track is calculated in a similar way to the `snpA` calculations:

`ra` = pa*(1-erate) + (1-pa)*erate
`rb` = (1-pa)*(1-erate) + pa*erate 

The overall `ra` and `rb` values are the product across all SNPs considered. These are then normalised to sum to 1.0.



<a name="unzipping"></a>
### Haplotig unzipping

Once phasing has been performed, the final step is to "unzip" the contigs into their haplotigs. SAMPhaser does this fairly crudely, and so it is important to do additional read mapping and error correction after unzipping. Only phase blocks with at least `minsnp=X` SNPs will be unzipped. Regions that are not unzipped are output as "Collapsed" (`C`) blocks.

The unzip process works by:

1. Identifying all the reads associated with a given phase block track. (The `*.haprid.tdt` file.)

2. Filtering haplotigs are filtered on the basis of depth of coverage. For each haplotig, the mean depth of coverage (e.g. summed read lengths divided by haplotig block length) is calculated. If this depth of coverage is less than `minhap=X` (default=5) then that haplotig is discarded.

3. Following haplotig filtering, there is the option to drop "orphan" haplotigs that are no longer phased into two haplotypes. This is controlled by `halfhap=T/F`. By default (True), "half haplotigs" are retained. **NOTE:** Setting `halfhap=F` might leave missing coverage in the final output. It is primarily used when output needs matched A & B haplotigs for assessment of phasing.

4. Extension. If a haplotig starts or terminates near the very end of the contig, it will be extended. This is controlled by the `endmargin=X` setting. The default is 10, i.e. a block will be extended by up to 10 nt to reach the start or end of the phased contig.

5. "Collapsed" (C) blocks are established based on the reads that do not fall into phased blocks. All the RIDs associated with a given locus are identified and reduced to the set that are not assigned to a phased track. Unassigned reads are then checked against phase blocks and any reads wholly contained within a phased haplotig (`Start` to `End`) are removed. **Note:** The full length is used for this, not just the phased region. This is to prevent short reads that fall outside the phased SNPs adding lots of small "C" regions that wholly overlap with phased A/B tracks.

6. Unzipped sequences are generated and output to `*.haplotigs.fas`. The reference sequence is used for the basis of the these sequences. Each allele for each SNP is considered and the dominant allele for a given haplotype is selected. If it has at least `unzipcut=X` (and `absunzipcut=X` if `unzipcut` < 1) `QN` reads supporting it then it is assigned to the haplotig. The final haplotig sequence is the original contig sequence with any assigned non-reference alleles substituted in at the appropriate positions. Deletions (`-` alleles) are cut out of the sequence and so it may end up shorter than the original contig. Insertions (`+N` and longer deletions `-N`) are not currently handled and are ignored; for this reason, it is important to re-map reads and correct the final haplotig sequences, e.g. with Quiver.

7. Finally, the subsets of reads assigned to haplotypes are used to generate new haplotype-specific read depth data. If `rgraphics=T` these will also be plotted for each locus. Because not all reads can be confidently assigned to a phased track, the summed depth of `A`+`B` tracks may not reach the total read depth for a given part of the original contig. Note also that positions in these plots are relative to the original contig, not the haplotig sequences, which may be shortened due to deletions.

Due to structural variants and rearrangements, it is possible for one haplotig to completely lose coverage for large chunks. These will be split into multiple haplotigs, controlled by the `splitzero=X` setting, where `X` is the minumum length (bp) of a region of 0X coverage required to split a haplotig. Setting `splitzero=-1` will disable splitting.

<a name="report"></a>
### HTML report generation

SNPs are plotted in the frequency and colour of the track that is different from the Reference sequence. If both alleles are non-Reference, the combined frequency is plotted in the colour of Track C.

_Details of the report will be added here._
