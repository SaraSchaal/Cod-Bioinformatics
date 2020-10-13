## Information about testing of --disable_adapter_trimming

There are summary output files for two samples (i.e., Pop3_17304 and Pop3_17325).
Each sample was run with --disable_adapter_trimming on with files that just have the sample name (e.g., files with Pop3_17304) and off that include "test" in the file name, (e.g., files with Pop3_17304_test). 
There does seem to be differences in these, but nothing like the one I looked at when I practiced. I wish I would have saved that practice output, but it basically had dropped most of the reads when it was turned on. Based on the output from Stacks I would think I would need to add this flag because the adapter seems to be trimmed already.