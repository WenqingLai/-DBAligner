/*
 * Configure.h
 *
 *  Created on: Mar 9, 2014
 *      Author: laiwenqi
 */

#ifndef CONFIGURE_H_
#define CONFIGURE_H_
#include<string>
using namespace std;
struct FileFormat{
	static int fa;
	static int fq;
	static int unknown;
};
class Configure{
public:
	static string config_file;
	static string dir;
	static string refs_file;
    static string reads_file;
    static string left_reads,right_reads;
	static string align_file;
	static int seed_size;
	static int max_refs;
	static int max_read_len;
	static float load_factor;
	static int seed_max_poses;
	static float global_edit_ratio,global_gap_ratio,local_edit_ratio;
	static int max_local_results;
	static int local_max_gap,local_min_edit;
	static int global_max_gap,global_max_edit;
	static float vote_cutoff;
	static bool report_all;
	static int ref_fmt;
	static int reads_fmt;
	static bool is_pe;
	static int insert_estimate_cutoff;
	static int backward_chimera_pos,forward_chimera_pos;
	static float deviation_factor;
	static string read1,read2;
	static string rev_read1,rev_read2;
	static float median_factor;
	static int low_bound,high_bound;
	static bool bound_designated;
public:
	static bool initialize();
	static bool initialize(const string & file);
	static void printSingleConfiguration();
	static void printPEConfiguration();
private:
	static int findReferenceFormat(const string & ref_file);
	static int findReadsFormat(const string & reads_file);
};


#endif /* CONFIGURE_H_ */
