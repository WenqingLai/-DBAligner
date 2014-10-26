/*
 * Graph.h
 *
 *  Created on: Mar 10, 2014
 *      Author: laiwenqi
 */

#ifndef GRAPH_H_
#define GRAPH_H_
#include "HashTable.h"
#include "Reference.h"
#include "Configure.h"
#include "Sam.h"
using namespace std;
struct Vertex{
	int content;
	unsigned char in_count,out_count;
	int in_edges[4],out_edges[4];
	static unsigned char set_flag;
	static unsigned char clear_flag;
};
struct SearchNode{
	int prev_pos;
	int vertex_loc;
	short read_pos;
	unsigned char action;
	unsigned char global_edit;
	unsigned char local_gap;
	unsigned char global_gap;
	unsigned char local_edit;
	static int no_prev;
	static unsigned char read_gap;
	static unsigned char path_gap;
	static unsigned char mismatch;
	static unsigned char match;
	static unsigned char start;
	static unsigned char terminal;
};
enum SearchKind{
	forward,middle,backward
};
class SearchInterval{
public:
	SearchKind type;
	int path_start,path_end;
	short read_start,read_end;
	short max_edit;
public:
	static int makeIntervals(int lookupShared[],bool marks[],int marks_len,SearchInterval intervals[]);
	static int makeLimitIntervals(int lookupShared[],bool marks[],int marks_len,SearchInterval intervals[]);
	static int makeBackwardIntervals(int lookupShared[],bool marks[],int marks_len,SearchInterval intervals[]);
};
struct NaiveAlignmentFormat{
	int ref_pos;
	short edit;
	short ref_read_len;
	string cigar;
};
class Graph {
private:
	string read;
	Vertex * vertices;
	int vertices_size,vertices_capacity;
	HashTable * pLookupString;
	int * lookupShared;
	bool * marks,* marks_mem;
	int marks_len;
	Reference * pRef_seq;
	const char * ref_str;
	int ref_size;
	SearchNode * searchNodes;
	int searchNodes_size,searchNodes_capacity;
	string part_read_buf,read_buf,part_cigar_buf,cigar_buf;
	SearchInterval * intervals;
	int intervals_size;
	int * valid_poses;
	int valid_poses_size;
	int * result_poses;
	int result_size;
	int * locations;
	int locations_capacity,locations_size;
	int approximate_start;
	NaiveAlignmentFormat * temp_alignments;
	int temp_alignments_capacity,temp_alignments_size;
	int * temp_alignments_poses;
	int * chimera_poses;
	int * pLookupInterval;
	int * best_poses;
	int best_poses_size;
public:
	AlignmentFormat alignment;
public:
	Graph(Reference * pRef_seq,HashTable * pHash_table);
	virtual ~Graph();
	int searchGraph(const string & read,bool temp_marks[],int len);
	int limitSearch(const string & read,int ref_start);
	int limitSearch2(const string & read,int ref_start);
	bool resolveAlignment();
	bool resolveLimitAlignment();
	bool resolveLimitAlignment(NaiveAlignmentFormat & n_alignment);
	int resolveAlignments(NaiveAlignmentFormat * alignments);
	int resolveBestAlignments(NaiveAlignmentFormat * alignments);
	int findStoneAlignment();
	int forwardChimera(NaiveAlignmentFormat * n_alignments,const string & read);
	int backwardChimera(NaiveAlignmentFormat * n_alignments,const string & read);
	void makeMarks(const string & read,bool temp_marks[],int len);
private:
	void constructGraph();
	bool findSharedVertex();
	void recifySharedVertex();
	bool findLimitSharedVertex();
	void backwardSearch(int active_start,const SearchInterval & interval);
	void middleSearch(int active_start,const SearchInterval & interval);
	void forwardSearch(int active_start,const SearchInterval & interval);
	void backwardLimitSearch(int active_start,const SearchInterval & interval);
	void middleLimitSearch(int active_start,const SearchInterval & interval);
	void forwardLimitSearch(int active_start,const SearchInterval & interval);
	void backwardLimitSearch2(int active_start,const SearchInterval & interval);
	void middleLimitSearch2(int active_start,const SearchInterval & interval);
	void forwardLimitSearch2(int active_start,const SearchInterval & interval);
	void sortTempAlignments(int a,int b);
	void sortAlignments(int a,int b);
	inline void clearLocations();
	void findBestPositions();
	void forwardChimeraSearch(int active_start,const SearchInterval & interval);
	void backwardChimeraSearch(int active_start,const SearchInterval & interval);
	void middleChimeraSearch(int active_start,const SearchInterval & interval);
	void backward_forwardChimeraSearch(int active_start,const SearchInterval & interval);
	void backward_backwardChimeraSearch(int active_start,const SearchInterval & interval);
	void backward_middleChimeraSearch(int active_start,const SearchInterval & interval);
	int resolveForwardChimeraAlignment(NaiveAlignmentFormat * n_alignments);
	int resolveBackwardChimeraAlignment(NaiveAlignmentFormat * n_alignments);
	void printIntervals();
	void printMarks();
	void printLocations();
	void printTempAlignments();
};
#endif /* GRAPH_H_ */
