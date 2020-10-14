## Email exchange with Jon Puritz about stacks, fastp and bwa

# My email:

	Hey Jon,

	Here are the things you requested in the meeting last week.

	Github account for dDocent version from oyster genome: SaraSchaal ( https://github.com/SaraSchaal )
	Example of how the reads look pre and post Stacks ( https://github.com/SaraSchaal/Cod-Bioinformatics/blob/master/notebook/20201007_Pre%26PostStacksExample.md )
	Example of FastP out files with –disable_adapter_trimming on and off for two samples: ( https://github.com/SaraSchaal/Cod-Bioinformatics/tree/master/src/FastP/20201007_testDisableTrimming )
	Bed file with repeats for Cod Genome. This doesn’t seem to exists on the NCBI website. There is a GFF file that I can potentially search for repeat elements in, but there doesn’t seem to be a separate file with repeats. ( https://www.ncbi.nlm.nih.gov/genome/?term=txid8049[orgn] )

	The adapter trimming is confusing to me. I couldn’t recreate what happened when I was practicing FastP (banging my head for not documenting those outputs that I got). However, based on the output from stacks, it seems like the adapters are trimmed. We get different outputs though with the flag added or not.  So I am curious to hear your thoughts on that. Thanks
	 
	The last thing is I am getting set up to run BWA but a couple of options I was wondering if you had any advice on:

	-A matching score (defaults to 1)

	-O gap penalty (defaults to 6)

	-B mismatch penalty (defaults to 4)
	 
	I have looked up what these are and how they work, but I’m not sure how to choose the right values for each one. I’ll keep digging online, but if you have any advice on those that would be appreciated!

	All the best,
	Sara

# Jon's Response:

	Hi Sara,

	Let me just respond to all of your items in turn:

	You’re now added. You should get an email for you to accept the invitation.
	
	The output post process_radtags looks good to me. You may want to consider adding the “rescue” flag -r to have the program rescue any barcodes that are off by one bp. This would increase your reads post demultiplexing (but only about 1% or so). I’d also advise you to look at the output log file. I have found with custom adapters that sometimes the barcode is off by 1 bp in the 5’ or 3’ direction. For example, the barcode ACTTGA might appear as CTTGA or TACTTGA. This should again be a small percentage of reads, but it’s worth checking on.
	
	So, process_radtags is simply removing your barcode and only if it’s perfectly at the beginning of your read (which it should be). However, removing adapters in fastp is about removing additional sequencing adapter that appears in your read. This happens when your insert length is shorter than your read length, so that end of the sequencing read goes past the insert and starts to go into the downstream (or upstream) sequencing adapter. Process_radtags will not catch these, at least not nearly as well as fastp.
	
	For the repeat bed file, take a look on the FTP site (https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Gadus_morhua/latest_assembly_versions/GCF_902167405.1_gadMor3.0/). The file GCF_902167405.1_gadMor3.0_rm.out.gz is going to have your identified repeats via the NCBI annotation pipeline.
	
	For bwa, I typically adjust the gap penalty to 5 and the mismatch penalty to 3. I find that this consistently gives me better results with a lot of marine species. However, it’s worth doing a couple of test runs with 3–5 samples comparing that to the default. This is an old exercise that’s specific to RADseq, but it has all of the concepts and ways to evaluate different mapping parameters (https://github.com/jpuritz/Winter.School2018/blob/master/Exercises/Day1/Mapping%20Exercise.md).

	My general advice for choosing parameter values is to choose them by testing, not based on defaults or (too much) advice online.

	Jon