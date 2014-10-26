/*
 * Sam.h
 *
 *  Created on: Mar 19, 2014
 *      Author: laiwenqi
 */

#ifndef SAM_H_
#define SAM_H_
#include <fstream>
using namespace std;
struct AlignmentFormat{
	int ref_index;
	int start_pos;
	string cigar;
	int edit;
};
struct PEAlignmentFormat{
	int ref_indices[2];
	int start_poses[2];
	string cigar[2];
	int edits[2];
	int insert_size;
};
class Sam {
private:
	ofstream output;
	static int unmapped;
	static int mapped_neg;
	static int pos_secondary;
	static int neg_secondary;
	static int pe_positive_mapped_f;
	static int pe_positive_mapped_s;
	static int pe_negative_mapped_f;
	static int pe_negative_mapped_s;
	static int pe_disc_ff_f;
	static int pe_disc_fr_f;
	static int pe_disc_rf_f;
	static int pe_disc_rr_f;
	static int pe_disc_ff_s;
	static int pe_disc_fr_s;
	static int pe_disc_rf_s;
	static int pe_disc_rr_s;
public:
	Sam();
	virtual ~Sam();
	bool writeOpen(const char * fn);
	void writeHeader(int ref_lens[],int size);
	void writeUnmappedAlignment(int read_num,const string & read);
	void writePosMappedAlignment(int read_num,int edit,const string & cigar,const string & read,int * prefs,int * poses,int size);
	void writeNegMappedAlignment(int read_num,int edit,const string & cigar,const string & read,int * prefs,int * poses,int size);
	void writePosPerfectAlignment(int read_num,const string & read,int * prefs,int * poses,int size);
	void writeNegPerfectAlignment(int read_num,const string & read,int * prefs,int * poses,int size);
	void writePosMappedAlignment(int read_num,const string & read,const AlignmentFormat & alignment);
	void writeNegMappedAlignment(int read_num,const string & read,const AlignmentFormat & alignment);
	void writePosSecMappedAlignment(int read_num,const AlignmentFormat & alignment);
	void writeNegSecMappedAlignment(int read_num,const AlignmentFormat & alignment);
	void writePosPEAlignment(int read_num,const PEAlignmentFormat & pe_alignment);
	void writeNegPEAlignment(int read_num,const PEAlignmentFormat & pe_alignment);
	void writePEDiscordance(int read_num,const PEAlignmentFormat & pe_alignment,bool r1_positive,bool r2_positive);
	void writePEUnmappedAlignment(int read_num,const PEAlignmentFormat & pe_alignment,bool r1_positive,bool r2_positive);
	void writePEUnmappedAlignment(int read_num);
	void writeClose();
private:
	string compressCigar(const string & raw_cigar);
	string compressCigar(const string * raw_cigar);
};

#endif /* SAM_H_ */
