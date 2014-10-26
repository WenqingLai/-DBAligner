/*
 * HashTable.h
 *
 *  Created on: Apr 21, 2014
 *      Author: laiwenqi
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_
#include "Reference.h"
#include "Configure.h"
struct Element{
	int ref_pos;
	int vertex_loc;
	int next;
};
struct FinalElement{
	int ref_pos;
	int next;
};
class HashTable {
private:
	const char * ref_str;
	int ref_size;
	int * pTable;
	int table_size;
	int * repeats;
	Element * buckets;
	int buckets_size;
	int buckets_capacity;
	FinalElement * presult;
	int result_capacity;
public:
	HashTable(Reference * pRef_seq);
	virtual ~HashTable();
	int verifyResult(const string & str,int poses[]);
	int lookupString(const string & str,int poses[],int max_size);
	bool tryToAdd(int ref_pos,int vertex_loc);
	int get(int ref_pos);
	int get(const string & str);
	bool containsKey(const string & str);
private:
	int getRepeatHead(const string & str);
	int hashCode(int begin);
	int hashCode(const string & str);
	inline bool refEqual(int ref_pos1,int ref_pos2);
	inline bool refEqual2(const string & str,int ref_pos);
};

#endif /* HASHTABLE_H_ */
