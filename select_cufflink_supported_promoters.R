library(GenomicFeatures)
library(Biostrings)
library(rtracklayer)

# The transcripts.gtf from cufflink output
cufflink.gtf = 'transcripts.gtf'
# The tmap file from cuffcompare output
cuffcmp.tmap = 'cuffcmp.transcripts.gtf.tmap'
# The fasta of the genome
genome.fasta = '../NO_DESC_Thalassiosira_pseudonana.ASM14940v1.17.dna.toplevel.fa'
# The bases before and after TSS
bp.before.TSS = 100
bp.after.TSS = 0

# FPKM filter
min.fpkm = 10
# Class filter
class = c('=', 'c', 'j')

# Extract the transcript_id from the attributes of each row of a GTF
get.transcript.id = function(x) {
  x = unlist(strsplit(x,'; '))
  x = x[grep('transcript_id',x)]
  x = gsub('transcript_id \"','',x)
  x = gsub('\"','',x)
  unlist(x)
}

# Load the transcripts.gtf from cufflink output as table
gtf = read.table(file=cufflink.gtf,sep='\t',comment.char='',quote='')
colnames(gtf) = c('seqid','source','type','start','end','score','strand','phase','attribute')
# Take out unstranded transcripts from the GTF
gtf = gtf[gtf$strand %in% c('+','-'),]

# Load the tmap file from cuffcompare output as a table
tmap = read.table(file=cuffcmp.tmap,sep='\t',comment.char='',quote='',head=T)

# Filter on FPKM
sel = tmap[tmap$FPKM >= min.fpkm,]
# Filter on classes
sel = sel[sel$class_code %in% class,]

# Extract transcript_id from the GTF working on attribute
attribute = as.character(gtf$attribute)
attribute = strsplit(attribute,'; ')
transcript.id = lapply(attribute,get.transcript.id)

# Prepare a data frame from the GTF containing only transcript_id, chromosome, start, end, strand information
gtf$transcript.id = transcript.id
range.tab = gtf[gtf$type == 'transcript', c('transcript.id','seqid','start','end','strand')]

# Merge the cufflink GTF and the cuffcompare table for the selected transcripts and write a table
res = merge(sel,range.tab,by.x='cuff_id',by.y='transcript.id',sort=F)
cufflink.tab.sel = paste('SELECTED_TAB',cufflink.gtf,sep='_')

# Calculate promoters coordinates
pos = res$strand == '+'
promoter.start = as.vector(ifelse(pos,res$start-bp.before.TSS,res$end-bp.after.TSS))
promoter.end = as.vector(ifelse(pos,res$start+bp.after.TSS,res$end+bp.before.TSS))
promoter.strand = ifelse(pos,'+','-')
a=cbind(res,promoter.start,promoter.end,promoter.strand)

write.table(res,file=cufflink.tab.sel,sep='\t',quote=F,row.names=F,col.names=T)

# Save filtered GTF for visualization (colum 1:9 because we have added the 10th with the transcript.id)
cufflink.gtf.sel = paste('SELECTED_GTF',cufflink.gtf,sep='_')
gtf = gtf[gtf$transcript.id %in% sel$cuff_id,c(1:9)]
write.table(gtf,file=cufflink.gtf.sel,sep='\t',quote=F,row.names=F,col.names=F)

# Load the genome sequences to extract the promoter sequences
genome = readDNAStringSet(file=genome.fasta)







// a = subseq(genome['1'],start=1,end=10)
// complement(reverse(a))

//Priority Code Description
//= Complete match of intron chain
//c Contained 
//j Potentially novel isoform (fragment): at least one splice junction is shared with a reference transcript 
//e Single exon transfrag overlapping a reference exon and at least 10 bp of a reference intron, indicating a possible pre-mRNA fragment. 
//i A transfrag falling entirely within a reference intron 
//o Generic exonic overlap with a reference transcript 
//p Possible polymerase run-on fragment (within 2Kbases of a reference transcript) 
//r Repeat. Currently determined by looking at the soft-masked reference sequence and applied to transcripts where at least 50% of the bases are lower case 
//u Unknown, intergenic transcript 
//x Exonic overlap with reference on the opposite strand 
//s An intron of the transfrag overlaps a reference intron on the opposite strand (likely due to read mapping errors) 
//. (.tracking file only, indicates multiple classifications) 
//
//We will consider: '=', 'c', 'j'
//
//transdb = makeTranscriptDbFromGFF('transcripts.gtf',
//format='gtf',
//exonRankAttributeName='exon_number',
//dataSource='cuffcompare',
//species='T pseudonana')
//
//saveDb(transdb,file='TP.cuffcmp.accurate.sqlite')
