# Refining the polyadenylation site (PAS) of transcripts and expression quantification based on outputs of IsoQuant

## Requirements
Data preparing:  
| Data | Note |
|:-------:|:-------:|
| reference_genome.fasta | Reference genome |
| gene_annotation.gtf | Reference gene annotation |
| bam file | Aligned and sorted BAM file |

Softwares:
| Tool | Source |
|:-------:|:-------:|
| IsoQuant  | https://github.com/ablab/IsoQuant?tab=readme-ov-file#examples |
| gtfToGenePred | https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred |
| genePredToBed | https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/genePredToBed |
| trmap | https://github.com/gpertea/gffcompare|

## Running
1. Run isoquant
```
isoquant --reference /PATH/TO/reference_genome.fasta --genedb /PATH/TO/gene_annotation.gtf \
--complete_genedb --bam sample_name.bam --data_type nanopore --prefix sample_name -o isoquant_output_dir
```
- ```sample_name.bam```: aligned and sorted BAM file
-  ```isoquant_output_dir```: directory containing outputs of isoquant
2. Convert extended_annotation.gtf to bed
```
cd isoquant_output_dir/sample_name

# convert gtf to genePred
gtfToGenePred sample_name.extended_annotation.gtf sample_name.extended_annotation.genePred

# convert genePred to bed12
genePredToBed sample_name.extended_annotation.genePred sample_name.extended_annotation.bed12

# mark gene id and gene name in bed file
awk '
BEGIN { OFS = FS = "\t" }
ARGIND==1 {
    if ($3 == "transcript") {
        n = split($9, parts, ";")   
        for (i = 1; i <= n; i++) {
            split(parts[i], subparts, "\"")
            gsub(/^[ \t]+|[ \t]+$/, "", subparts[1])
            gsub(/^[ \t]+|[ \t]+$/, "", subparts[2])
            if (subparts[1] == "transcript_id") {
                transcript_id = subparts[2]
            }
            gsub(/^[ \t]+|[ \t]+$/, "", subparts[1])
            if (subparts[1] == "gene_id") {
                gene_id = subparts[2]
            }
            if (subparts[1] == "gene_name") {
                gene_name = subparts[2]
            }
        }
        dict[transcript_id] = transcript_id"|"gene_id"|"gene_name
    }        
}
ARGIND==2 {
    if ($4 in dict) {
        $4 = dict[$4]
        print $0
    }
}
' sample_name.extended_annotation.gtf \
sample_name.extended_annotation.bed12 \
> sample_name.extended_annotation.bed
```
3. Refine PAS
```
python add_pas_from_isoquant.py \
--annotation_bed sample_name.extended_annotation.bed \
--reads_assignments sample_name.read_assignments.tsv.gz \
--corrected_reads sample_name.corrected_reads.bed.gz \
--output_path pas_refined_dir
```
- ```sample_name.read_assignments.tsv.gz```：isoquant output TSV file with read to isoform assignments 
- ```sample_name.corrected_reads.bed.gz```: isoquant output BED file with corrected read alignments
-  ```pas_refined_dir```: directory containing PAS refined bed file
4. Extract full-length reads and quantify gene/isoform expression
```
gunzip sample_name.corrected_reads.bed.gz
# assigning reads with trmap
trmap pas_refined_dir/pas_refined.bed sample_name.corrected_reads.bed -o trmap_output_dir/sample_name.read_to_isoform.txt

# extract full-length reads and quantify gene/isoform expression
python reads_assignment.py \
--trmap trmap_output_dir/sample_name.read_to_isoform.txt \
--corrected_reads_bed sample_name.corrected_reads.bed \
--prefix sample_name \
-o read_assignment_output_dir
```
-  ```trmap_output_dir```: directory containing outputs of trmap
-  ```read_assignment_output_dir```: directory containing outputs of reads_assignment.py
