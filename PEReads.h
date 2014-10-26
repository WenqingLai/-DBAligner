/*
 * PEReads.h
 *
 *  Created on: May 8, 2014
 *      Author: laiwenqi
 */

#ifndef PEREADS_H_
#define PEREADS_H_
#include "Graph.h"
#include "SubRead.h"
#include "DNA.h"
using namespace std;
enum AlignedType{
	pe_mapped,pe_discordant,pe_chimeric,pe_unmapped
};
class PEReads {
private:
	int * inserts;
	int inserts_size,inserts_capacity,half_inserts_size;
	Reference * pRef_seq;
	Graph * pGraph;
	SubRead * pSub_read;
	bool * pmarks,* nmarks;
	NaiveAlignmentFormat * left_alignments, * right_alignments, * alignments;
	int left_size,right_size,alignments_capacity;
	int * alignments_poses;
public:
	PEAlignmentFormat pe_alignment;
	bool left_positive,right_positive;
	int average,st_deviation;
	int low_bound,high_bound;
public:
	PEReads(Reference * pRef_seq,Graph * pGraph, SubRead * pSub_read);
	virtual ~PEReads();
	bool estimateInsert();
	void caculateInsert();
	AlignedType alignPERead();
private:
	void calculateAverage(int start,int end);
	void calculateStDeviation(int start,int end);
	int alignSingleRead(const string & read,NaiveAlignmentFormat * alignments,bool & positive);
	int alignSingleStoneRead(const string & read);
	void sortAlignments(int start,int end);
	int lookupLowerBound(NaiveAlignmentFormat * alignments,int size,int ref_pos);
	int lookupUpperBound(NaiveAlignmentFormat * alignments,int size,int ref_pos);
	bool findPositiveBestPair();
	bool findNegativeBestPair();
	void printInserts();
	void printNaiveAlignments(NaiveAlignmentFormat * alignments,int size);
	void printPEAlignment(const PEAlignmentFormat & alignment);
};

#endif /* PEREADS_H_ */
