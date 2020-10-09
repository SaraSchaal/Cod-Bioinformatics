##### Notes for quality control #####
#### Trimmomatic ####
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

#### Fast P ####
### Notes on fastp for quality control (Jon P uses in new version of dDocent)
	fastp --in1 $1.F.fq.gz --in2 $1.R.fq.gz --out1 $1.R1.fq.gz --out2 $1.R2.fq.gz --cut_front --cut_tail --cut_window_size 5 --cut_mean_quality 15 --correction $TW -q 15 -u 50 -j $1.json -h $1.html --detect_adapter_for_pe &> $1.trim.log 

### run on example file
	fastp --in1 Stacks_Out/Pop3_18304.1.fq.gz --in2 Stacks_Out/Pop3_18304.2.fq.gz --out1 FastP_Out/Pop3_18304.1.R1.fq.gz --out2 FastP_Out/Pop3_18304.1.R2.fq.gz -q 15 -u 50 --trim_front1 1 --cut_front --cut_tail --disable_adapter_trimming --cut_window_size 5 --cut_mean_quality 15 -j Pop3_18304.1.json -h FastP_Out/Pop3_18304.1.html &> FastP_Out/Pop3_18304.1.trim.log

# ABOVE LINE OF CODE WORKED on personal computer: 20200720

### loop for practice on local drive
file = fastP_practice.sh
<!-- for i in Stacks_Out/*head.1.fq.gz 
do
	sampleID=${i:11:10};
	echo $sampleID;
	fastp --in1 Stacks_Out/$sampleID.head.1.fq.gz --in2 Stacks_Out/$sampleID.head.2.fq.gz --out1 FastP_Out/$sampleID.head.R1.fq.gz --out2 FastP_Out/$sampleID.head.R2.fq.gz -q 15 -u 50 --trim_front1 1 --cut_front --cut_tail --disable_adapter_trimming --cut_window_size 5 --cut_mean_quality 15 -j $sampleID.fp.json -h FastP_Out/$sampleID.fp.html &> FastP_Out/$sampleID.fp.trim.log
done 
 -->

# ABOVE LINE OF CODE WORKED on personal computer: 20200819

### Order of files to run on cluster
1) create submission files 
2) loop through submission files and submit to cluster