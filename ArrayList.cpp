/*
 * ArrayList.cpp
 *
 *  Created on: Apr 30, 2014
 *      Author: laiwenqi
 */

#include "ArrayList.h"
#include "Configure.h"
#include <iostream>

ArrayList::ArrayList(int capacity) {
	list_capacity=capacity;
	pArrayList=new ArrayListElement[list_capacity];
	ref_starts=new int[Configure::max_local_results];
}

ArrayList::~ArrayList() {
	// TODO Auto-generated destructor stub
}
bool ArrayList::addRefPositions(int poses[],int size,int count,int offset){//size >=1
	int i,j,ref_start,prev_i;
	ArrayListElement element;
	element.next_pos=-1;
	if(list_size==0){
		element.ref_start=poses[0]-offset;
		element.count=count;
		pArrayList[0]=element;
		list_size=1;
		for(i=1;i!=size;i++){
			pArrayList[list_size-1].next_pos=list_size;
			element.ref_start=poses[i]-offset;
			element.count=count;
			pArrayList[list_size++]=element;
			if(list_size==list_capacity)
				return false;
		}
		head=0;
		return true;
	}
	i=head;j=0;prev_i=-1;
	while(i!=-1&&j!=size){
		ref_start=poses[j]-offset;
		if(pArrayList[i].ref_start<ref_start){
			prev_i=i;
			i=pArrayList[i].next_pos;
		}
		else if(pArrayList[i].ref_start>ref_start){
			element.ref_start=ref_start;
			element.count=count;
			element.next_pos=i;
			pArrayList[list_size]=element;
			if(prev_i!=-1)
				pArrayList[prev_i].next_pos=list_size;
			else
				head=list_size;
			list_size++;
			if(list_size==list_capacity)
				return false;
			j++;
		}
		else{
			pArrayList[i].count+=count;
			prev_i=i;
			i=pArrayList[i].next_pos;
			j++;
		}
	}
	for(;j!=size;j++){
		pArrayList[prev_i].next_pos=list_size;
		element.ref_start=poses[j]-offset;
		element.count=count;
		prev_i=list_size;
		pArrayList[list_size++]=element;
		if(list_size==list_capacity)
			return false;
	}
	pArrayList[prev_i].next_pos=-1;
	return true;
}
void ArrayList::getAllPossibleRefStarts(int cutoff){
	int prev_i,i,prev_value,vote;
	prev_value=pArrayList[head].ref_start;
	vote=pArrayList[head].count;
	prev_i=head;
	i=pArrayList[head].next_pos;
	ref_starts_size=0;
	while(i!=-1){
		if(pArrayList[i].ref_start-prev_value<=Configure::global_max_gap)
			vote+=pArrayList[i].count;
		else{
			if(vote>=cutoff){
				ref_starts[ref_starts_size++]=prev_value;
				if(ref_starts_size==Configure::max_local_results)
					break;
			}
			prev_value=pArrayList[i].ref_start;
			vote=pArrayList[i].count;
		}
		prev_i=i;
		i=pArrayList[i].next_pos;
	}
	if(ref_starts_size!=Configure::max_local_results&&pArrayList[prev_i].ref_start-prev_value<=Configure::global_max_gap&&vote>=cutoff){
		ref_starts[ref_starts_size++]=prev_value;
	}
}
int ArrayList::getAllPossibleRefStarts(int cutoff,int candidate_poses[]){
	int prev_i,i,prev_value,vote,poses_size;
	prev_value=pArrayList[head].ref_start;
	vote=pArrayList[head].count;
	prev_i=head;
	i=pArrayList[head].next_pos;
	poses_size=0;
	while(i!=-1){
		if(pArrayList[i].ref_start-prev_value<=Configure::global_max_gap)
			vote+=pArrayList[i].count;
		else{
			if(vote>=cutoff){
				candidate_poses[poses_size++]=prev_value;
				if(poses_size==Configure::max_local_results)
					break;
			}
			prev_value=pArrayList[i].ref_start;
			vote=pArrayList[i].count;
		}
		prev_i=i;
		i=pArrayList[i].next_pos;
	}
	if(poses_size!=Configure::max_local_results&&pArrayList[prev_i].ref_start-prev_value<=Configure::global_max_gap&&vote>=cutoff){
		candidate_poses[poses_size++]=prev_value;
	}
	return poses_size;
}
int ArrayList::getMostPossibleRefStart2(){
	int i,max_vote,prev_value,vote,max_vote_pos;
	max_vote=-1;
	prev_value=pArrayList[head].ref_start;
	vote=pArrayList[head].count;
	i=pArrayList[head].next_pos;
	while(i!=-1){
		if(pArrayList[i].ref_start-prev_value<=Configure::global_max_gap)
			vote+=pArrayList[i].count;
		else{
			if(vote>max_vote){
				max_vote=vote;
				max_vote_pos=prev_value;
			}
			prev_value=pArrayList[i].ref_start;
			vote=pArrayList[i].count;
		}
		i=pArrayList[i].next_pos;
	}
	if(vote>max_vote){
		max_vote_pos=prev_value;
		max_vote=vote;
	}
	return max_vote_pos;
}
void ArrayList::getMostPossibleRefStart(){
	int i,max_vote,prev_value,vote,max_vote_pos;
	max_vote=-1;
	prev_value=pArrayList[head].ref_start;
	vote=pArrayList[head].count;
	i=pArrayList[head].next_pos;
	while(i!=-1){
		if(pArrayList[i].ref_start-prev_value<=Configure::global_max_gap)
			vote+=pArrayList[i].count;
		else{
			if(vote>max_vote){
				max_vote=vote;
				max_vote_pos=prev_value;
			}
			prev_value=pArrayList[i].ref_start;
			vote=pArrayList[i].count;
		}
		i=pArrayList[i].next_pos;
	}
	if(vote>max_vote){
		max_vote_pos=prev_value;
		max_vote=vote;
	}
	ref_starts_size=0;
	ref_starts[ref_starts_size++]=max_vote_pos;
}
void ArrayList::printArrayList(){
	int i;
	i=head;
	while(i!=-1){
		cout<<pArrayList[i].ref_start<<","<<pArrayList[i].count<<" ";
		i=pArrayList[i].next_pos;
	}
	cout<<endl;
}
