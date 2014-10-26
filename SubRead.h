/*
 * SubRead.h
 *
 *  Created on: May 1, 2014
 *      Author: laiwenqi
 */

#ifndef SUBREAD_H_
#define SUBREAD_H_

#include "Reference.h"
#include "Sam.h"
#include "HashTable.h"
#include "ArrayList.h"
#include "Graph.h"
class SubRead{
private:
	int ref_size;
	const char * ref_str;
	HashTable * pTable;
	ArrayList * pRef_starts;
	Graph * pGraph;
	Reference * pRef_seq;
	int ref_starts_capacity;
	int * seed_poses, * left_seed_poses,* right_seed_poses;
	int seed_poses_size,left_seed_size,right_seed_size;
	AlignmentFormat * temp_alignments;
	int temp_alignments_size;
	int * temp_alignments_poses;
	bool * marks;
	int marks_len;
public:
	AlignmentFormat ** alignments;
	int alignments_size;
	AlignmentFormat alignment;
	PEAlignmentFormat pe_alignment;
public:
	SubRead(Reference * pRef_seq,HashTable * pTable,Graph * pGraph);
	virtual ~SubRead();
	bool alignRead(const string & read,bool marks[],int marks_len);
	int alignRead(const string & read,bool marks[],int marks_len,NaiveAlignmentFormat * alignments);
	int alignStoneRead(const string & read,bool marks[],int marks_len);
	bool alignPERead(int low_bound,int high_bound);
private:
	void printMarks(bool marks[],int len);
	void sortAlignments(int start,int end);
	int findSinglePositions(const string & read,int seed_poses[]);
	int lookupLowerBound(int poses[],int size,int goal);
	int lookupUpperBound(int poses[],int size,int goal);
	bool findPositiveBestPair(int low_bound,int high_bound);
	bool findNegativeBestPair(int low_bound,int high_bound);
	void recifyMarks();
};

#endif /* SUBREAD_H_ */
