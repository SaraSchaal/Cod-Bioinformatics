### Notes on Trimmomatic for quality control 
Leading - cut bases off the start of a read, if below a threshold quality
Headcrop - cut the specified number of bases from the start of the read
Slidingwindow - window size is number of bases to average across and required quality specifies the average quality required

# Should we cut the first base to not worry about the A overhang. I wasn't finding any papers that discuss this. Can reach out to Jon if you don't think this is a good thing to do without knowing more?
# What quality should we trim?
	Q=30, fails 1 in 1000, Q=20, fails 1 in 100. 
	so should we do Q > 20?
	See what Jon does in dDocent (sliding window of 5 basepairs), and look at recent cod papers. 
	A overhang just chop the first base with headcrop. 

# line of code to manipulate for paired end reads
	java -jar <path to trimmomatic.jar> PE [-threads <threads] [-phred33 | -phred64] [-trimlog <logFile>] <input 1> <input 2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <step 1>
# NOT USING TRIMMOMATIC JON'S NEW PIPELINE USES FASTP ALL INPUTS ARE THE SAME 

### Notes on fastp for quality control (Jon P uses in new version of dDocent)
	fastp --in1 $1.F.fq.gz --in2 $1.R.fq.gz --out1 $1.R1.fq.gz --out2 $1.R2.fq.gz --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction $TW -q 15 -u 50 -j $1.json -h $1.html --detect_adapter_for_pe &> $1.trim.log 

### run on example file
	fastp --in1 Stacks_Out/Pop3_18304.1.fq.gz --in2 Stacks_Out/Pop3_18304.2.fq.gz --out1 FastP_Out/Pop3_18304.1.R1.fq.gz --out2 FastP_Out/Pop3_18304.1.R2.fq.gz -q 15 -u 50 --trim_front1 1 --cut_front --cut_tail --disable_adapter_trimming --cut_window_size 5 --cut_mean_quality 15 -j Pop3_18304.1.json -h Pop3_18304.1.html &> Pop3_18304.1.trim.log

### Examples for stepping through files 
#Trimmomatic - NOT USING ANYMORE
<!-- for infile in *_1.fastq.gz
 do
   base=$(basename ${infile} _1.fastq.gz)
   trimmomatic PE ${infile} ${base}_2.fastq.gz \
                ${base}_1.trim.fastq.gz ${base}_1un.trim.fastq.gz \
                ${base}_2.trim.fastq.gz ${base}_2un.trim.fastq.gz \
                SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15 
 done -->


<!-- for i in `ls -1 *1.fq.gz | sed 's/\_R1.fastq//'`; 
do 
	echo trimmomatic PE -phred33 $i\_1.fq.gz $i\_2.fq.gz $i\_R1_paired.fq.gz $i\_R1_unpaired.fq.gz $i\_R2_paired.fq.gz $i\_R2_unpaired.fq.gz ILLUMINACLIP:contams_forward_rev.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> cmd_file; 
done -->

#FastP
for i in `ls -1 *1.fq.gz | sed 's/\_R1.fastq//'`;

	fastp --in1 $1.fq.gz --in2 $2.fq.gz --out1 $1.R1.fq.gz --out2 $1.R2.fq.gz --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 -q 15 -u 50 --disable_adapter_trimming --trim_front1 1 -j $1.json -h $1.html &> $1.trim.log 

done
