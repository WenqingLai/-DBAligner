-DBAligner
==========

Novel short reads aligner via De Bruijn graph approach

This aligner has a higher mapping ratio and accuracy within similar time consumption when tested on bacterial genome, compared with Bowtie, Bowtie2 and BWA.
The memory footprint is large,about 64X reference size,and it's not proper for large and complex genome like the whole human genome.

Just download the package and complile it in Linux environment.Makefile is hand-coded and it's very simple to make. 
To work in Windows, just modify the "Main.cpp" file and delete "#include <sys/resource.h> #include <sys/time.h>" and the corresponding statements.

The behavior of DBAligner can be manipulated by editing the "runtime_parameters" file and DBAligner can read in parameters from this file rather than from the terminal.
To launch DBAligner, just type "DBAligner" in the terminal if the file "runtime_parameters" loactes in the same directory, or type "DBAligner runtime_parameters" intead if not.

Runtime parameters are divided mainly to three categories: IO, single reads and paired reads.

IO:
##data_dir="./"
Input and output directory
##ref_file="reference.fa"
If reference has many sequences, it should be merged in one fa file
##reads_file="left_reads.fq right_reads.fq" 
Or reads_file="reads.fq" in single reads alignment. Reads file can also be fa format
##aln_result="DBAligner.sam"
Aln_result can only be SAM format.

Single reads:
##seed_size=18
Seed_size is tuned not too high and not too low for accuracy vs speed
##global_edits_ratio=0.15
If read length is 100bp, the allowed edits is 15bp.
##global_gap_bases_ratio=0.08
If read length is 100bp, the allowed gap bases is 8bp. 
Caution,this is different from other aligners. We control gap bases rather than gap numbers.
##local_edits_ratio=0.32
One read can partially aligned in many places. The unmapped parts form extension intervals.
We control edits in each extension interval.
##local_min_edits=3
If the extension interval is 2bp,the local_edits_ratio may derive 0.64bp edits.
Maximum_Allowed_Edits=Max{local_min_edits,local_edits_ratio*extension_interval_len}
##local_max_gap_bases=2
Likewise,we control the gap bases in each extension interval.
For efficiency,we don't let the allowed local gap bases grow along with the extension interval length as ##local_edits_ratio has done.
##max_local_results=200
This controls the local temporary alignments in each extension interval.
##report_all=0
If the value is 1, all alignments are reported in SAM file.
##max_refs=200
The maximum allowed sequences in reference.For instance, the human genome only need 23 sequences or more if the sequences are not continuous.
##max_read_len=500
We can align reads with various lengths ranging from seed size to max_read_len.

Paired reads exclusive:
##vote_cutoff=0.2
When either end of pe can be mapped.DBAligner realign the pe by aligning the two ends simultanously.
This parameter controls candidate positions selection. 
For instance, one end of pe can be sheared into 100 seeds and only 50 can be mapped.
Of the 50 seeds,only 9 indicates position x, 41 indicates position y.
Position x is filtered out because there's only 0.18 support ratio.
We don't extend alignment from position x to save much time.
##lest_stone_pes=100
This parameter deals with insert size estimation.
The insert size can be derived if both ends of pe can by mapped once and only once.
We call the above alignment stone alignment.
So,if 200 pes are mapped, but only 99 stone alignments, we need to align more pes to estimate insert size.
##median_factor=0.5
This parameter also deals with insert size estimation.
There are 100 stone alignments and 100 insert sizes are derived.
We choose the median to signify the estimated insert size rather than the average value.
If median_factor=0.45, the estimation is inclined to a smaller value.
##deviation_factor=4
Insert size obey normal distribution. Deviation_factor can be 3 to include more than 0.90 cases.
##designate_insert_scope=0
DBAligner can estimate insert size and standard deviation if this parameter is set to 0
Otherwise, you must give the insert size scope by the following two parameters
##insert_low_bound=0
The smallest allowed insert size.
##insert_high_bound=600
The largest allowed insert size.





