## 20210223 - email from Shobana about FreeBayes

Hi Sara,

Thanks for opening the ticket and reaching out to us! Do let us know how your meeting with the collaborator goes, but in the meantime, here are a few suggestions you may consider:

1) Split your run into different chromosomes and you can later concatenate the results together. Freebayes does have an option to run on specific chr as below:
freebayes -f ref.fa -r chrQ aln.bam >var.vcf

2) freebayes-parallel - https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel - this will enable you to run freebayes in parallel on for e.g., 100000bp chunks of the reference
https://github.com/freebayes/freebayes/blob/master/scripts/freebayes-parallel

3) If you have any sample data you can share with us, we can try and help you write a shell script to concatenate and post-process the results. 

4) Lastly, you may always split your jobs between the short and lotterhos partitions to enable maximum utilization of resources.

Do let us know if your collaborator had additional suggestions too.

Thanks,
Shobana
