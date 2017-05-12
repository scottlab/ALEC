# ALEC 
Amplicon Long-read Error Correction (ALEC) was developed to correct sequencing and alignment errors [substitutionsandinsertion/deletions (indels)] generated by targeted amplicon sequencing with the PacBio RS platform. This script has been developed usingPacBio single molecule real time (SMRT) full-gene sequencing of the CYP2D6 gene (5.0 kb)according to the P6-C4 Pacific Biosciences protocol. ALEC was further tested using a 9.2kb amplicon and long-read PacBio SMRT sequencing of the RYR2 gene. ALEC may have utility withother long-read sequencing platforms (e.g., Oxford Nanopore), and/or other sequencing chemistries; however, optimizing the error correction parameters may be necessary to achieve ideal sequence correction and output for genes and sequencing platforms other than those noted above.

Four types of sequencing errors can be corrected after applying ALEC: 1) random substitutions; 2) random indels; 3) indels within homopolymers; and 4) indels near sequence variants.

 
## System requirements
* Python 2.7.10 or above   
* Does nott support Python 3 
* pysam
* numpy

## Usage
The ALEC script takes a fasta file as a reference file and a SAM/BAM file as the raw data alignment file, and automatically generates a corrected fasta file as output in the working directory. The usage of the script is as below:

Python ALEC.py -r reference.fasta -i input.bam/sam [arguments]

### Arguments Table

|Argument|Type|Default|Description|
|:-----|:---:|:---:|:---|
|--input <br/>   -i	| string	| NA	| Required. Input file, SAM or BAM. |
|-- reference<br/> -r|string|	NA	|Required. Reference file, FASTA file, index it following Samtools manual| 
|--targetRegion<br/> -t|	string|	NA	|Required. Target region interval. Example: 1:300000-400000|
|--lengthFilter<br/> -lf	|float|	0.001|	Optional. Reads shorter than lengthFilter * length(targetRegion)<br/> will be excluded in this correction process.|
|--downsample<br/>	-ds|	float|	1.0|	Optional. Fraction of down sampling. |
|--deletion<br/> -del|	float|	0.0	|Required. Deletion error frequency threshold (per base) to trigger correction.|
|--insert<br/> -ins|	float|	0.0|	Required. Insert error frequency threshold (per base) to trigger correction.|
|--mismatch<br/> -mis|	float	|0.0|	Required. Substitution error frequency threshold (per base) to trigger correction.|
|--del_homo_p<br/>-del_hp |	float	|0.0|	Required. Deletion Homopolymer Penalty.|
|--ins_homo_p<br/> -ins_hp |	float	|0.0|	Optional. Insert Homopolymer Penalty.|
|--platform<br/> -x|string|NA|Optional. Use preset arguments for correction.|
|--help <br/> -h| | |show help message|

## Note
1. The script only takes one single sequence as reference each time. Please use the same reference file as used in alignment. 
2. We do not have a preference regarding the available sequence alignment tools; however, we used alignment files generated with BWA-MEM (0.7.12) to develop the ALEC script (see reference below for details). Files from other alignment tool could lead to unexpected performance. 

## Reference
1.  A manuscript detailing and evaluating the functionality of ALEC is currently under review.  
2.  For an application example of ALEC, see reference: 
  Qiao W and Yang Y, et al. Hum Mutat. 2016 Mar;37(3):315-23. doi: 10.1002/humu.22936. Epub 2015 Dec 18.
  
## Contact  
* [SCOTT Lab at Mount Sinai](http://stuartscottlab.org/)
