/*
 * ReadMap.cpp
 *
 *  Created on: Apr 30, 2014
 *      Author: laiwenqi
 */

#include "SubRead.h"
#include  "DNA.h"
#include <cstdlib>
#include <algorithm>
#include <iostream>
SubRead::SubRead(Reference * pRef_seq,HashTable * pTable, Graph * pGraph) {
	ref_size=pRef_seq->ref_str.size();
	ref_str=pRef_seq->ref_str.c_str();
	this->pTable=pTable;
	this->pRef_seq=pRef_seq;
	ref_starts_capacity=1<<12;
	pRef_starts=new ArrayList(ref_starts_capacity);
	this->pGraph=pGraph;
	temp_alignments=new AlignmentFormat[Configure::max_local_results];
	temp_alignments_poses=new int[Configure::max_local_results];
	alignments=new AlignmentFormat*[Configure::max_local_results];
	left_seed_poses=new int[Configure::max_local_results];
	right_seed_poses=new int[Configure::max_local_results];
	seed_poses=left_seed_poses;
	marks=new bool[Configure::max_read_len];
}
SubRead::~SubRead() {
}
void SubRead::printMarks(bool marks[],int marks_size){
	int i;
	for(i=0;i!=marks_size;i++)
		if(marks[i])
			cout<<"1";
		else
			cout<<"0";
	cout<<endl;
}
bool SubRead::alignRead(const string & read,bool marks[],int marks_size){
	int i,count,start,temp_size,vote_cutoff;
	int * ref_starts,valid_size,vote_kmers;
	bool in_gap;
	in_gap=true;
	vote_kmers=0;
	pRef_starts->list_size=0;
	for(i=0;i!=marks_size;i++){
		if(marks[i]){
			if(in_gap){
				count=1;
				in_gap=false;
			}
			else
				count++;
		}
		else if(!in_gap){
			in_gap=true;
			start=i-count;
			temp_size=count+Configure::seed_size-1;
			seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
			if(seed_poses_size==0){
				if(count>1){
					start=i-count+(rand()%count)/2;
					temp_size=(count+1)/2+Configure::seed_size-1;
					seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
					if(seed_poses_size==0)
						continue;
					else
						vote_kmers+=(count+1)/2;
				}
				else
					continue;
			}
			else
				vote_kmers+=count;
			if(!pRef_starts->addRefPositions(seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start))
				break;
		}
	}
	if(!in_gap){
		in_gap=true;
		start=i-count;
		temp_size=count+Configure::seed_size-1;
		seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
		if(seed_poses_size==0){
			if(count>1){
				start=i-count+(rand()%count)/2;
				temp_size=(count+1)/2+Configure::seed_size-1;
				seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
				if(seed_poses_size!=0)
					vote_kmers+=(count+1)/2;
			}
		}
		else
			vote_kmers+=count;
		if(seed_poses_size!=0)
			pRef_starts->addRefPositions(seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start);
	}
	if(pRef_starts->list_size==0)
		return false;
	vote_cutoff=Configure::vote_cutoff*vote_kmers;
	pRef_starts->getAllPossibleRefStarts(vote_cutoff);
	if(pRef_starts->ref_starts_size==0){
		pRef_starts->getMostPossibleRefStart();
		if(pRef_starts->list_size==0)
			return false;
	}
	temp_size=pRef_starts->ref_starts_size;
	ref_starts=pRef_starts->ref_starts;
	alignment.edit=0x7fffffff;
	for(i=0;i!=temp_size;i++){
		valid_size=pGraph->limitSearch2(read,ref_starts[i]);
		if(valid_size==0)
			continue;
		if(!pGraph->resolveLimitAlignment())
			continue;
		if(alignment.edit>pGraph->alignment.edit)
			alignment=pGraph->alignment;
	}
	if(alignment.edit==0x7fffffff)
		return false;
	return true;
}
int SubRead::alignRead(const string & read,bool marks[],int marks_size,NaiveAlignmentFormat * alignments){
	int i,count,start,temp_size,vote_cutoff,alignments_size;
	int * ref_starts,valid_size,vote_kmers;
	bool in_gap;
	NaiveAlignmentFormat n_alignment;
	in_gap=true;
	vote_kmers=0;
	pRef_starts->list_size=0;
	for(i=0;i!=marks_size;i++){
		if(marks[i]){
			if(in_gap){
				count=1;
				in_gap=false;
			}
			else
				count++;
		}
		else if(!in_gap){
			in_gap=true;
			start=i-count;
			temp_size=count+Configure::seed_size-1;
			seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
			if(seed_poses_size==0){
				if(count>1){
					start=i-count+(rand()%count)/2;
					temp_size=(count+1)/2+Configure::seed_size-1;
					seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
					if(seed_poses_size==0)
						continue;
					else
						vote_kmers+=(count+1)/2;
				}
				else
					continue;
			}
			else
				vote_kmers+=count;
			if(!pRef_starts->addRefPositions(seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start))
				break;
		}
	}
	if(!in_gap){
		in_gap=true;
		start=i-count;
		temp_size=count+Configure::seed_size-1;
		seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
		if(seed_poses_size==0){
			if(count>1){
				start=i-count+(rand()%count)/2;
				temp_size=(count+1)/2+Configure::seed_size-1;
				seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
				if(seed_poses_size!=0)
					vote_kmers+=(count+1)/2;
			}
		}
		else
			vote_kmers+=count;
		if(seed_poses_size!=0)
			pRef_starts->addRefPositions(seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start);
	}
	if(pRef_starts->list_size==0)
		return 0;
	vote_cutoff=Configure::vote_cutoff*vote_kmers;
	pRef_starts->getAllPossibleRefStarts(vote_cutoff);
	if(pRef_starts->ref_starts_size==0){
		pRef_starts->getMostPossibleRefStart();
		if(pRef_starts->list_size==0)
			return 0;
	}
	temp_size=pRef_starts->ref_starts_size;
	ref_starts=pRef_starts->ref_starts;
	alignments_size=0;
	for(i=0;i!=temp_size;i++){
		valid_size=pGraph->limitSearch2(read,ref_starts[i]);
		if(valid_size==0)
			continue;
		if(pGraph->resolveLimitAlignment(n_alignment))
			alignments[alignments_size++]=n_alignment;
		if(alignments_size==Configure::max_local_results)
			return alignments_size;
	}
	return alignments_size;
}
int SubRead::alignStoneRead(const string & read,bool marks[],int marks_size){
	int i,count,start,temp_size,vote_cutoff;
	int * ref_starts,valid_pos_size,vote_kmers;
	int ref_pos,temp_ref_pos;
	bool in_gap,mapped;
	in_gap=true;
	vote_kmers=0;
	pRef_starts->list_size=0;
	for(i=0;i!=marks_size;i++){
		if(marks[i]){
			if(in_gap){
				count=1;
				in_gap=false;
			}
			else
				count++;
		}
		else if(!in_gap){
			in_gap=true;
			start=i-count;
			temp_size=count+Configure::seed_size-1;
			seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
			if(seed_poses_size==0){
				if(count>1){
					start=i-count+(rand()%count)/2;
					temp_size=(count+1)/2+Configure::seed_size-1;
					seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
					if(seed_poses_size==0)
						continue;
					else
						vote_kmers+=(count+1)/2;
				}
				else
					continue;
			}
			else
				vote_kmers+=count;
			if(!pRef_starts->addRefPositions(seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start))
				break;
		}
	}
	if(!in_gap){
		in_gap=true;
		start=i-count;
		temp_size=count+Configure::seed_size-1;
		seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
		if(seed_poses_size==0){
			if(count>1){
				start=i-count+(rand()%count)/2;
				temp_size=(count+1)/2+Configure::seed_size-1;
				seed_poses_size=pTable->lookupString(read.substr(start,temp_size),seed_poses,Configure::seed_max_poses);
				if(seed_poses_size!=0)
					vote_kmers+=(count+1)/2;
			}
		}
		else
			vote_kmers+=count;
		if(seed_poses_size!=0)
			pRef_starts->addRefPositions(seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start);
	}
	if(pRef_starts->list_size==0)
		return -1;
	vote_cutoff=Configure::vote_cutoff*vote_kmers;
	pRef_starts->getAllPossibleRefStarts(vote_cutoff);
	if(pRef_starts->ref_starts_size==0){
		pRef_starts->getMostPossibleRefStart();
		if(pRef_starts->list_size==0)
			return -1;
	}
	temp_size=pRef_starts->ref_starts_size;
	ref_starts=pRef_starts->ref_starts;
	temp_alignments_size=0;
	mapped=false;
	for(i=0;i!=temp_size;i++){
		valid_pos_size=pGraph->limitSearch2(read,ref_starts[i]);
		if(valid_pos_size==0)
			continue;
		temp_ref_pos=pGraph->findStoneAlignment();
		if(temp_ref_pos==-2)
			return -2;
		else if(temp_ref_pos==-1)
			continue;
		if(!mapped){
			mapped=true;
			ref_pos=temp_ref_pos;
		}
		else if(abs(ref_pos-temp_ref_pos)>Configure::global_max_gap){
			return -2;
		}
	}
	if(mapped)
		return ref_pos;
	else
		return -1;
}
void SubRead::sortAlignments(int start,int end){
	int i,j,pivot;
	pivot=temp_alignments_poses[start];
	i=start;j=end-1;
	while(i<j){
		while(i<j&&temp_alignments[temp_alignments_poses[j]].edit>temp_alignments[pivot].edit)
			j--;
		if(i==j)
			break;
		temp_alignments_poses[i]=temp_alignments_poses[j];
		i++;
		while(i<j&&temp_alignments[temp_alignments_poses[i]].edit<temp_alignments[pivot].edit)
			i++;
		if(i==j)
			break;
		temp_alignments_poses[j]=temp_alignments_poses[i];
		j--;
	}
	temp_alignments_poses[i]=pivot;
	if(start<i)
	   sortAlignments(start,i);
	if(i+1<end)
	   sortAlignments(i+1,end);
}
bool SubRead::alignPERead(int low_bound,int high_bound){
	left_seed_size=findSinglePositions(Configure::read1,left_seed_poses);
	right_seed_size=findSinglePositions(Configure::rev_read2,right_seed_poses);
	if(left_seed_size&&right_seed_size&&findPositiveBestPair(low_bound,high_bound))
		return true;
	left_seed_size=findSinglePositions(Configure::rev_read1,left_seed_poses);
	right_seed_size=findSinglePositions(Configure::read2,right_seed_poses);
	if(left_seed_size&&right_seed_size&&findNegativeBestPair(low_bound,high_bound))
		return true;
	return false;
}
int SubRead::findSinglePositions(const string & read,int temp_seed_poses[]){
	int i,vote_kmers,count,start,temp_size,vote_cutoff;
	bool in_gap;
	marks_len=read.length()-Configure::seed_size+1;
	for(i=0;i!=marks_len;i++){
		if(pTable->containsKey(read.substr(i,Configure::seed_size)))
			marks[i]=true;
		else
			marks[i]=false;
	}
	recifyMarks();
	in_gap=true;
	vote_kmers=0;
	pRef_starts->list_size=0;
	for(i=0;i!=marks_len;i++){
		if(marks[i]){
			if(in_gap){
				count=1;
				in_gap=false;
			}
			else
				count++;
		}
		else if(!in_gap){
			in_gap=true;
			start=i-count;
			temp_size=count+Configure::seed_size-1;
			seed_poses_size=pTable->lookupString(read.substr(start,temp_size),temp_seed_poses,Configure::seed_max_poses);
			if(seed_poses_size==0){
				if(count>1){
					start=i-count+(rand()%count)/2;
					temp_size=(count+1)/2+Configure::seed_size-1;
					seed_poses_size=pTable->lookupString(read.substr(start,temp_size),temp_seed_poses,Configure::seed_max_poses);
					if(seed_poses_size==0)
						continue;
					else
						vote_kmers+=(count+1)/2;
				}
				else
					continue;
			}
			else
				vote_kmers+=count;
			if(!pRef_starts->addRefPositions(temp_seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start))
				break;
		}
	}
	if(!in_gap){
		in_gap=true;
		start=i-count;
		temp_size=count+Configure::seed_size-1;
		seed_poses_size=pTable->lookupString(read.substr(start,temp_size),temp_seed_poses,Configure::seed_max_poses);
		if(seed_poses_size==0){
			if(count>1){
				start=i-count+(rand()%count)/2;
				temp_size=(count+1)/2+Configure::seed_size-1;
				seed_poses_size=pTable->lookupString(read.substr(start,temp_size),temp_seed_poses,Configure::seed_max_poses);
				if(seed_poses_size!=0)
					vote_kmers+=(count+1)/2;
			}
		}
		else
			vote_kmers+=count;
		if(seed_poses_size!=0)
			pRef_starts->addRefPositions(temp_seed_poses,seed_poses_size,temp_size-Configure::seed_size+1,start);
	}
	if(pRef_starts->list_size==0)
		return 0;
	vote_cutoff=Configure::vote_cutoff*vote_kmers;
	temp_size=pRef_starts->getAllPossibleRefStarts(vote_cutoff,temp_seed_poses);
	if(temp_size==0){
		temp_seed_poses[0]=pRef_starts->getMostPossibleRefStart2();
		return 1;
	}
	return temp_size;
}
int SubRead::lookupLowerBound(int poses[],int size,int goal){
	int start,end,mid;
	start=-1;end=size;
	while(start<end-1){
		mid=(start+end)/2;
		if(goal<=poses[mid])
			end=mid;
		else
			start=mid;
	}
	return start;
}
int SubRead::lookupUpperBound(int poses[],int size,int goal){
	int start,end,mid;
	start=-1;end=size;
	while(start<end-1){
		mid=(start+end)/2;
		if(goal<poses[mid])
			end=mid;
		else
			start=mid;
	}
	return end;
}
bool SubRead::findPositiveBestPair(int low_bound,int high_bound){
	int i,j,start,end,valid_size;
	int read1_len,read2_len;
	NaiveAlignmentFormat n_alignment1,n_alignment2;
	read1_len=Configure::read1.length();
	read2_len=Configure::read2.length();
	if(left_seed_size<right_seed_size){
		if(left_seed_size>1)
			sort(left_seed_poses,left_seed_poses+left_seed_size);
		for(i=0;i!=right_seed_size;i++){
			start=lookupLowerBound(left_seed_poses,left_seed_size,right_seed_poses[i]-high_bound+read2_len);
			end=lookupUpperBound(left_seed_poses,left_seed_size,right_seed_poses[i]-low_bound+read2_len);
			if(start+1>=end)
				continue;
			for(j=start+1;j!=end;j++){
				valid_size=pGraph->limitSearch2(Configure::read1,left_seed_poses[j]);
				if(valid_size&&pGraph->resolveLimitAlignment(n_alignment1))
					break;
			}
			if(j==end)
				continue;
			valid_size=pGraph->limitSearch2(Configure::rev_read2,right_seed_poses[i]);
			if(valid_size&&pGraph->resolveLimitAlignment(n_alignment2)){
				pe_alignment.start_poses[0]=n_alignment1.ref_pos;
				pe_alignment.start_poses[1]=n_alignment2.ref_pos;
				pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len);
				if(pe_alignment.ref_indices[0]!=pe_alignment.ref_indices[1])
					continue;
				pe_alignment.cigar[0]=n_alignment1.cigar;
				pe_alignment.cigar[1]=n_alignment2.cigar;
				pe_alignment.edits[0]=n_alignment1.edit;
				pe_alignment.edits[1]=n_alignment2.edit;
				pe_alignment.insert_size=n_alignment2.ref_pos-n_alignment1.ref_pos+n_alignment2.ref_read_len;
				break;
			}
		}//for
		if(i==right_seed_size)
			return false;
		else
			return true;
	}//if
	if(right_seed_size>1)
		sort(right_seed_poses,right_seed_poses+right_seed_size);
	for(i=0;i!=left_seed_size;i++){
		start=lookupLowerBound(right_seed_poses,right_seed_size,left_seed_poses[i]+low_bound-read2_len);
		end=lookupUpperBound(right_seed_poses,right_seed_size,left_seed_poses[i]+high_bound-read2_len);
		if(start+1>=end)
			continue;
		for(j=start+1;j!=end;j++){
			valid_size=pGraph->limitSearch2(Configure::rev_read2,right_seed_poses[j]);
			if(valid_size&&pGraph->resolveLimitAlignment(n_alignment2))
				break;
		}
		if(j==end)
			continue;
		valid_size=pGraph->limitSearch2(Configure::read1,left_seed_poses[i]);
		if(valid_size&&pGraph->resolveLimitAlignment(n_alignment1)){
			pe_alignment.start_poses[0]=n_alignment1.ref_pos;
			pe_alignment.start_poses[1]=n_alignment2.ref_pos;
			pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len);
			if(pe_alignment.ref_indices[0]!=pe_alignment.ref_indices[1])
					continue;
			pe_alignment.cigar[0]=n_alignment1.cigar;
			pe_alignment.cigar[1]=n_alignment2.cigar;
			pe_alignment.edits[0]=n_alignment1.edit;
			pe_alignment.edits[1]=n_alignment2.edit;
			pe_alignment.insert_size=n_alignment2.ref_pos-n_alignment1.ref_pos+n_alignment2.ref_read_len;
			break;
		}
	}//for
	if(i==left_seed_size)
		return false;
	else
		return true;
}
bool SubRead::findNegativeBestPair(int low_bound,int high_bound){
	int i,j,start,end,valid_size;
	int read1_len,read2_len;
	NaiveAlignmentFormat n_alignment1,n_alignment2;
	read1_len=Configure::rev_read1.length();
	read2_len=Configure::read2.length();
	if(left_seed_size<right_seed_size){
		if(left_seed_size>1)
			sort(left_seed_poses,left_seed_poses+left_seed_size);
		for(i=0;i!=right_seed_size;i++){
			start=lookupLowerBound(left_seed_poses,left_seed_size,right_seed_poses[i]+low_bound+read1_len);
			end=lookupUpperBound(left_seed_poses,left_seed_size,right_seed_poses[i]+high_bound+read1_len);
			if(start+1>=end)
				continue;
			for(j=start+1;j!=end;j++){
				valid_size=pGraph->limitSearch2(Configure::rev_read1,left_seed_poses[j]);
				if(valid_size&&pGraph->resolveLimitAlignment(n_alignment1))
					break;
			}
			if(j==end)
				continue;
			valid_size=pGraph->limitSearch2(Configure::read2,right_seed_poses[i]);
			if(valid_size&&pGraph->resolveLimitAlignment(n_alignment2)){
				pe_alignment.start_poses[0]=n_alignment1.ref_pos;
				pe_alignment.start_poses[1]=n_alignment2.ref_pos;
				pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len);
				if(pe_alignment.ref_indices[0]!=pe_alignment.ref_indices[1])
					continue;
				pe_alignment.cigar[0]=n_alignment1.cigar;
				pe_alignment.cigar[1]=n_alignment2.cigar;
				pe_alignment.edits[0]=n_alignment1.edit;
				pe_alignment.edits[1]=n_alignment2.edit;
				pe_alignment.insert_size=n_alignment2.ref_pos-n_alignment1.ref_pos-n_alignment1.ref_read_len;
				break;
			}
		}//for
		if(i==right_seed_size)
			return false;
		else
			return true;
	}//if
	if(right_seed_size>1)
		sort(right_seed_poses,right_seed_poses+right_seed_size);
	for(i=0;i!=left_seed_size;i++){
		start=lookupLowerBound(right_seed_poses,right_seed_size,left_seed_poses[i]-high_bound+read1_len);
		end=lookupUpperBound(right_seed_poses,right_seed_size,left_seed_poses[i]-low_bound+read1_len);
		if(start+1>=end)
			continue;
		for(j=start+1;j!=end;j++){
			valid_size=pGraph->limitSearch2(Configure::read2,right_seed_poses[j]);
			if(valid_size&&pGraph->resolveLimitAlignment(n_alignment2))
				break;
		}
		if(j==end)
			continue;
		valid_size=pGraph->limitSearch2(Configure::rev_read1,left_seed_poses[i]);
		if(valid_size&&pGraph->resolveLimitAlignment(n_alignment1)){
			pe_alignment.start_poses[0]=n_alignment1.ref_pos;
			pe_alignment.start_poses[1]=n_alignment2.ref_pos;
			pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len);
			if(pe_alignment.ref_indices[0]!=pe_alignment.ref_indices[1])
					continue;
			pe_alignment.cigar[0]=n_alignment1.cigar;
			pe_alignment.cigar[1]=n_alignment2.cigar;
			pe_alignment.edits[0]=n_alignment1.edit;
			pe_alignment.edits[1]=n_alignment2.edit;
			pe_alignment.insert_size=n_alignment2.ref_pos-n_alignment1.ref_pos-n_alignment1.ref_read_len;
			break;
		}
	}//for
	if(i==left_seed_size)
		return false;
	else
		return true;
}
void SubRead::recifyMarks(){
	int i,island_count,gap_count;
	bool in_gap;
	gap_count=0;
	for(i=0;i!=marks_len;i++){
		if(marks[i])
			break;
		else
			gap_count++;
	}
	island_count=0;
	for(;i!=marks_len;i++){
		if(marks[i])
			island_count++;
		else
			break;
	}
	if(i==marks_len)
		return;
	if(gap_count==0){
		if(island_count>=2){
			marks[i-1]=false;
			gap_count=1;
		}
		else{
			gap_count=0;
		}
	}
	else{
		if(island_count>=3){
			marks[i-island_count]=false;
			marks[i-1]=false;
			gap_count=1;
		}
		else if(island_count==2){
			marks[i-1]=false;
			gap_count=1;
		}
		else{
			gap_count=0;
		}
	}
	in_gap=true;
	for(;i!=marks_len;i++){
		if(marks[i]){
			if(in_gap){
				in_gap=false;
				island_count=1;
			}
			else{
				island_count++;
			}
		}
		else{
			if(in_gap){
				gap_count++;
			}
			else{
				if(island_count>=3){
					marks[i-island_count]=false;
					marks[i-1]=false;
					if(gap_count+1<Configure::seed_size&&!marks[i-island_count-gap_count-2])
						marks[i-island_count-gap_count-1]=false;
					gap_count=2;
				}
				else if(island_count==2){
					if(gap_count+1<Configure::seed_size){
						marks[i-2]=false;
						marks[i-1]=false;
						gap_count+=3;
					}
					else if(gap_count<Configure::seed_size){
						marks[i-2]=false;
						gap_count=1;
					}
					else{
						marks[i-1]=false;
						gap_count=2;
					}
				}
				else if(island_count==1){
					if(gap_count<Configure::seed_size){
						marks[i-1]=false;
						gap_count+=2;
					}
					else
						gap_count=1;
				}
				in_gap=true;
			}
		}
	}
	if(!in_gap){
		if(island_count>=3){
			marks[i-island_count]=false;
			if(gap_count+1<Configure::seed_size&&!marks[i-island_count-gap_count-2])
					marks[i-island_count-gap_count-1]=false;
		}
		else if(island_count==2)
			marks[i-2]=false;
	}
}
