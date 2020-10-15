## Stacks with and without -r (rescue) flag
The “rescue” flag -r has the program rescue any barcodes that are off by one bp. This would increase your reads post demultiplexing (but only about 1% or so). I’d also advise you to look at the output log file. I have found with custom adapters that sometimes the barcode is off by 1 bp in the 5’ or 3’ direction. For example, the barcode ACTTGA might appear as CTTGA or TACTTGA. This should again be a small percentage of reads, but it’s worth checking on.

# with -r flag
	srun process_radtags -P -p GenomeFilesi59_i72/ -o ../Stacks_Out/rerun -b barcodeFilei59_i72.txt -r --inline_inline --disable_rad_check

Processing paired-end data.
Using Phred+33 encoding for quality scores.
Found 1 paired input file(s).
Searching for single and paired-end, inlined barcodes.
Loaded 4 barcodes (6bp / 6bp).
Will attempt to recover barcodes with at most 1 / 1 mismatches.
Processing file 1 of 1 [i5-9-i7-2_R1_001.fastq.gz]
  Reading data from:
  GenomeFilesi59_i72/i5-9-i7-2_R1_001.fastq.gz and
  GenomeFilesi59_i72/i5-9-i7-2_R2_001.fastq.gz
  Processing RAD-Tags...1M...2M...3M...4M...5M...6M...7M...8M...9M...10M...11M...12M...13M...14M...15M...16M...17M...18M...19M...20M...21M...22M...23M...24M...25M...26M...27M...28M...29M...30M...31M...32M...33M...34M...35M...36M...37M...38M...39M...40M...41M...42M...43M...44M...45M...46M...47M...48M...49M...50M...51M...52M...53M...54M...55M...56M...57M...58M...59M...60M...61M...62M...63M...64M...65M...66M...67M...68M...69M...70M...71M...72M...73M...74M...75M...76M...77M...78M...79M...80M...81M...82M...83M...84M...85M...86M...87M...88M...89M...90M...91M...92M...93M...94M...95M...96M...97M...98M...99M...100M...101M...102M...103M...104M...105M...106M...107M...108M...109M...110M...111M...112M...113M...114M...115M...116M...117M...118M...119M...120M...121M...122M...123M...124M...125M...126M...127M...128M...129M...130M...131M...132M...133M...134M...135M...136M...137M...138M...139M...140M...141M...142M...143M...144M...145M...146M...147M...148M...149M...150M...151M...152M...153M...154M...155M...156M...157M...158M...159M...160M...161M...162M...163M...164M...165M...166M...167M...168M...169M...170M...171M...172M...173M...174M...175M...176M...177M...178M...179M...180M...181M...182M...
  365761434 total reads; -13530898 ambiguous barcodes; -0 ambiguous RAD-Tags; +5721982 recovered; -0 low quality reads; 352230536 retained reads.
Closing files, flushing buffers...
Outputing details to log: '../Stacks_Out/rerun/process_radtags.GenomeFilesi59_i72.log'

365761434 total sequences
 13530898 barcode not found drops (3.7%)
        0 low quality read drops (0.0%)
        0 RAD cutsite not found drops (0.0%)
352230536 retained reads (96.3%)


# without -r flag
	srun process_radtags -P -p GenomeFilesi59_i72/ -o ../Stacks_Out -b barcodeFilei59_i72.txt --inline_inline --disable_rad_check

Processing paired-end data.
Using Phred+33 encoding for quality scores.
Found 1 paired input file(s).
Searching for single and paired-end, inlined barcodes.
Loaded 4 barcodes (6bp / 6bp).
Processing file 1 of 1 [i5-9-i7-2_R1_001.fastq.gz]
  Reading data from:
  GenomeFiles/i5-9-i7-2_R1_001.fastq.gz and
  GenomeFiles/i5-9-i7-2_R2_001.fastq.gz
  Processing RAD-Tags...1M...2M...3M...4M...5M...6M...7M...8M...9M...10M...11M...12M...13M...14M...15M...16M...17M...18M...19M...20M...21M...22M...23M...24M...25M...26M...27M...28M...29M...30M...31M...32M...33M...34M...35M...36M...37M...38M...39M...40M...41M...42M...43M...44M...45M...46M...47M...48M...49M...50M...51M...52M...53M...54M...55M...56M...57M...58M...59M...60M...61M...62M...63M...64M...65M...66M...67M...68M...69M...70M...71M...72M...73M...74M...75M...76M...77M...78M...79M...80M...81M...82M...83M...84M...85M...86M...87M...88M...89M...90M...91M...92M...93M...94M...95M...96M...97M...98M...99M...100M...101M...102M...103M...104M...105M...106M...107M...108M...109M...110M...111M...112M...113M...114M...115M...116M...117M...118M...119M...120M...121M...122M...123M...124M...125M...126M...127M...128M...129M...130M...131M...132M...133M...134M...135M...136M...137M...138M...139M...140M...141M...142M...143M...144M...145M...146M...147M...148M...149M...150M...151M...152M...153M...154M...155M...156M...157M...158M...159M...160M...161M...162M...163M...164M...165M...166M...167M...168M...169M...170M...171M...172M...173M...174M...175M...176M...177M...178M...179M...180M...181M...182M...
  365761434 total reads; -19252880 ambiguous barcodes; -0 ambiguous RAD-Tags; +0 recovered; -0 low quality reads; 346508554 retained reads.
Closing files, flushing buffers...
Outputing details to log: '../Stacks_Out/process_radtags.GenomeFiles.log'

365761434 total sequences
 19252880 barcode not found drops (5.3%)
        0 low quality read drops (0.0%)
        0 RAD cutsite not found drops (0.0%)
346508554 retained reads (94.7%)

# outcome: rerun all samples with rescue flag we retained 1.6% more reads by having this flag added
