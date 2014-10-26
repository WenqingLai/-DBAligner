/*
 * ArrayList.h
 *
 *  Created on: Apr 30, 2014
 *      Author: laiwenqi
 */

#ifndef ARRAYLIST_H_
#define ARRAYLIST_H_
struct ArrayListElement{
	int ref_start;
	int count;
	int next_pos;
};
class ArrayList {
private:
	ArrayListElement * pArrayList;
	int list_capacity,head;
public:
	int list_size;
	int * ref_starts;
	int ref_starts_size;
public:
	ArrayList(int capacity);
	virtual ~ArrayList();
	bool addRefPositions(int poses[],int size,int count,int offset);//poses must be in ascending order
	int getMostPossibleRefStart2();
	void getMostPossibleRefStart();
	void printArrayList();
	void getAllPossibleRefStarts(int cutoff);
	int getAllPossibleRefStarts(int cutoff,int candidate_poses[]);
};

#endif /* ARRAYLIST_H_ */
