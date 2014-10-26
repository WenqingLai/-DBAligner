/*
 * Reference.h
 *
 *  Created on: Apr 21, 2014
 *      Author: laiwenqi
 */

#ifndef REFERENCE_H_
#define REFERENCE_H_
#include<string>
#include "Configure.h"
#include "Sam.h"
using namespace std;
class Reference {
private:
	static char table[];
	static char mask;
public:
	string ref_str;
	int * refs_start;
	int refs_start_size;
	int total_N;
public:
	Reference();
	virtual ~Reference();
	bool purifyAndFillGap();
	int countRefIndexs(int * poses,const string & cigar,int * refs,int len,int read_len);
	int countRefIndexs(int * poses,int * refs,int len,int read_len);
	int countRefIndexs(int * poses,int delta,int * refs,int len);
	void countRefIndex(int & ref_index,int & start_pos,const string & cigar,int read_len);
	void countRefIndex(int & ref_index,int & start_pos,int read_len);
	void countRefIndex(AlignmentFormat & alignment,int read_len);
	bool countRefIndex(PEAlignmentFormat & pe_alignment,int read1_len, int read2_len);
private:
	void printRefsStart();
};

#endif /* REFERENCE_H_ */
