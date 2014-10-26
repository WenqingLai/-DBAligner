/*
 * PEReads.cpp
 *
 *  Created on: May 8, 2014
 *      Author: laiwenqi
 */

#include "PEReads.h"
#include<cmath>
#include<cstdlib>
#include<iostream>
#include<iomanip>

PEReads::PEReads(Reference * pRef_seq,Graph * pGraph,SubRead * pSub_read) {
	this->pRef_seq=pRef_seq;
	this->pGraph=pGraph;
	this->pSub_read=pSub_read;
	pmarks=new bool[Configure::max_read_len];
	nmarks=new bool[Configure::max_read_len];
	inserts_capacity=Configure::insert_estimate_cutoff;
	inserts_size=0;
	inserts=new int[inserts_capacity];
	alignments_capacity=Configure::max_local_results;
	left_alignments=new NaiveAlignmentFormat[alignments_capacity];
	right_alignments=new NaiveAlignmentFormat[alignments_capacity];
	alignments_poses=new int[alignments_capacity];
}

PEReads::~PEReads() {
}
bool PEReads::estimateInsert(){
	int i,read1_len,read2_len;
	int ref_pos1,ref_pos2;
	read1_len=Configure::read1.length();
	read2_len=Configure::read2.length();
	ref_pos1=alignSingleStoneRead(Configure::read1);
	if(ref_pos1<0)
		return false;
	ref_pos2=alignSingleStoneRead(Configure::read2);
	if(ref_pos2<0)
		return false;
	if(ref_pos1>ref_pos2)
		inserts[inserts_size++]=ref_pos1-ref_pos2+read1_len;
	else
		inserts[inserts_size++]=ref_pos2-ref_pos1+read2_len;
	if(inserts_size!=inserts_capacity)
		return false;
	half_inserts_size=inserts_size*Configure::median_factor;
	calculateAverage(0,inserts_size);
	for(i=0;i!=inserts_size;i++)
		inserts[i]=abs(inserts[i]-average);
	calculateStDeviation(0,inserts_size);
	low_bound=average-Configure::deviation_factor*st_deviation;
	if(low_bound<0)
		low_bound=0;
	high_bound=average+Configure::deviation_factor*st_deviation;
	return true;
}
void PEReads::caculateInsert(){
	int i;
	low_bound=high_bound=average=0;
	half_inserts_size=inserts_size*0.48;
	calculateAverage(0,inserts_size);
	for(i=0;i!=inserts_size;i++)
		inserts[i]=abs(inserts[i]-average);
	calculateStDeviation(0,inserts_size);
	low_bound=average-Configure::deviation_factor*st_deviation;
	if(low_bound<0)
		low_bound=0;
	high_bound=average+Configure::deviation_factor*st_deviation;
}
AlignedType PEReads::alignPERead(){
	int read1_len,read2_len,start_pos,ref_index;
	bool in_scope;
	NaiveAlignmentFormat alignment;
	PEAlignmentFormat temp_pe_alignment;
	left_size=alignSingleRead(Configure::read1,left_alignments,left_positive);
	right_size=alignSingleRead(Configure::read2,right_alignments,right_positive);
	temp_pe_alignment.ref_indices[0]=temp_pe_alignment.ref_indices[1]=-1;
	pe_alignment.ref_indices[0]=pe_alignment.ref_indices[1]=-1;
	if(left_size&&right_size){
		in_scope=false;
		if(left_positive&&!right_positive)
			in_scope=findPositiveBestPair();
		else if(!left_positive&right_positive)
			in_scope=findNegativeBestPair();
		if(in_scope)
			return pe_mapped;
	}
	if(pSub_read->alignPERead(low_bound,high_bound)){
		pe_alignment=pSub_read->pe_alignment;
		return pe_mapped;
	}
	if(!left_size&&!right_size){
		return pe_unmapped;
	}
	read1_len=Configure::read1.length();
	read2_len=Configure::read2.length();
	if(left_size&&!right_size){
		start_pos=left_alignments[0].ref_pos;
		pRef_seq->countRefIndex(ref_index,start_pos,left_alignments[0].ref_read_len);
		if(ref_index!=-1){
			temp_pe_alignment.ref_indices[0]=ref_index;
			temp_pe_alignment.start_poses[0]=start_pos;
			temp_pe_alignment.cigar[0]=left_alignments[0].cigar;
			temp_pe_alignment.edits[0]=left_alignments[0].edit;
		}
		Configure::global_max_edit=Configure::global_edit_ratio*read2_len;
		Configure::global_max_gap=Configure::global_gap_ratio*read2_len;
		if(left_positive){
			right_size=pGraph->forwardChimera(right_alignments,Configure::rev_read2);
			if(right_size==0){
				pe_alignment=temp_pe_alignment;
				return pe_unmapped;
			}
			in_scope=findPositiveBestPair();
		}
		else{
			right_size=pGraph->backwardChimera(right_alignments,Configure::read2);
			if(right_size==0){
				pe_alignment=temp_pe_alignment;
				return pe_unmapped;
			}
			in_scope=findNegativeBestPair();
		}
		if(in_scope)
			return pe_chimeric;
		else{
			pe_alignment=temp_pe_alignment;
			return pe_unmapped;
		}
	}
	if(right_size&&!left_size){
		start_pos=right_alignments[0].ref_pos;
		pRef_seq->countRefIndex(ref_index,start_pos,right_alignments[0].ref_read_len);
		if(ref_index!=-1){
			temp_pe_alignment.ref_indices[1]=ref_index;
			temp_pe_alignment.start_poses[1]=start_pos;
			temp_pe_alignment.cigar[1]=right_alignments[0].cigar;
			temp_pe_alignment.edits[1]=right_alignments[0].edit;
		}
		Configure::global_max_edit=Configure::global_edit_ratio*read1_len;
		Configure::global_max_gap=Configure::global_gap_ratio*read1_len;
		if(right_positive){
			left_size=pGraph->forwardChimera(left_alignments,Configure::rev_read1);
			if(left_size==0){
				pe_alignment=temp_pe_alignment;
				return pe_unmapped;
			}
			in_scope=findNegativeBestPair();
		}
		else{
			left_size=pGraph->backwardChimera(left_alignments,Configure::read1);
			if(left_size==0){
				pe_alignment=temp_pe_alignment;
				return pe_unmapped;
			}
			in_scope=findPositiveBestPair();
		}
		if(in_scope)
			return pe_chimeric;
		else{
			pe_alignment=temp_pe_alignment;
			return pe_unmapped;
		}
	}
	pe_alignment.start_poses[0]=left_alignments[0].ref_pos;
	pe_alignment.start_poses[1]=right_alignments[0].ref_pos;
	pe_alignment.cigar[0]=left_alignments[0].cigar;
	pe_alignment.cigar[1]=right_alignments[0].cigar;
	pe_alignment.edits[0]=left_alignments[0].edit;
	pe_alignment.edits[1]=right_alignments[0].edit;
	pe_alignment.insert_size=0;
	pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len);
	if(pe_alignment.ref_indices[0]==-1||pe_alignment.ref_indices[1]==-1)
		return pe_unmapped;
	return pe_discordant;
}
void PEReads::calculateAverage(int start,int end){
	int i,j,pivot;
	pivot=inserts[start];
	i=start;j=end-1;
	while(i<j){
		while(i<j&&inserts[j]>pivot)
			j--;
		if(i==j)
			break;
		inserts[i]=inserts[j];
		i++;
		while(i<j&&inserts[i]<pivot)
			i++;
		if(i==j)
			break;
		inserts[j]=inserts[i];
		j--;
	}
	inserts[i]=pivot;
	if(i==half_inserts_size){
		average=pivot;
		return;
	}
	else if(i<half_inserts_size){
		calculateAverage(i+1,end);
	}
	else{
		calculateAverage(start,i);
	}
}
void PEReads::calculateStDeviation(int start,int end){
	int i,j,pivot;
	i=rand()%(end-start)+start;
	pivot=inserts[i];
	inserts[i]=inserts[start];
	i=start;j=end-1;
	while(i<j){
		while(i<j&&inserts[j]>pivot)
			j--;
		if(i==j)
			break;
		inserts[i]=inserts[j];
		i++;
		while(i<j&&inserts[i]<pivot)
			i++;
		if(i==j)
			break;
		inserts[j]=inserts[i];
		j--;
	}
	inserts[i]=pivot;
	if(i==half_inserts_size){
		st_deviation=pivot;
		return;
	}
	else if(i<half_inserts_size){
		calculateStDeviation(i+1,end);
	}
	else{
		calculateStDeviation(start,i);
	}
}
int PEReads::alignSingleRead(const string & read,NaiveAlignmentFormat * alignments,bool & positive){
	int alignments_size,read_len,temp_marks_len;
	string rev_read;
	NaiveAlignmentFormat n_alignment;
	read_len=read.length();
	Configure::global_max_edit=read_len*Configure::global_edit_ratio;
	Configure::global_max_gap=read_len*Configure::global_gap_ratio;
	alignments_size=0;
	temp_marks_len=read_len-Configure::seed_size+1;
	pGraph->makeMarks(read,pmarks,temp_marks_len);
	if(pGraph->searchGraph(read,pmarks,temp_marks_len))
		alignments_size=pGraph->resolveBestAlignments(alignments);
	if(alignments_size){
		positive=true;
		return alignments_size;
	}
	rev_read=DNA::reverseComplete(read);
	pGraph->makeMarks(rev_read,nmarks,temp_marks_len);
	if(pGraph->searchGraph(rev_read,nmarks,temp_marks_len))
		alignments_size=pGraph->resolveBestAlignments(alignments);
	if(alignments_size){
		positive=false;
		return alignments_size;
	}
	return 0;
}
int PEReads::alignSingleStoneRead(const string & read){
	int valid_pos_size,read_len,prev_ref_pos,ref_pos,temp_marks_len;
	string rev_read;
	ref_pos=prev_ref_pos=-1;
	read_len=read.length();
	Configure::global_max_edit=read_len*Configure::global_edit_ratio;
	Configure::global_max_gap=read_len*Configure::global_gap_ratio;
	temp_marks_len=read_len-Configure::seed_size+1;
	pGraph->makeMarks(read,pmarks,temp_marks_len);
	valid_pos_size=pGraph->searchGraph(read,pmarks,temp_marks_len);
	if(valid_pos_size){
		ref_pos=pGraph->findStoneAlignment();
		if(ref_pos==-2)
			return ref_pos;
		prev_ref_pos=ref_pos;
	}
	rev_read=DNA::reverseComplete(read);
	pGraph->makeMarks(rev_read,nmarks,temp_marks_len);
	valid_pos_size=pGraph->searchGraph(rev_read,nmarks,temp_marks_len);
	if(valid_pos_size){
		ref_pos=pGraph->findStoneAlignment();
		if(ref_pos==-2)
			return ref_pos;
		if(ref_pos==-1)
			return prev_ref_pos;
		if(prev_ref_pos==-1)
			return ref_pos;
		if(abs(prev_ref_pos-ref_pos)>Configure::global_max_gap)
			return -2;
	}
	return ref_pos;
	//pGraph->copyMarks(nmarks,nmarks_len);
	//ref_pos=pSub_read->alignStoneRead(read,pmarks,pmarks_len);
	//if(ref_pos!=-1)
		//return ref_pos;
	//ref_pos=pSub_read->alignStoneRead(rev_read,nmarks,nmarks_len);
	//return ref_pos;
}
void PEReads::sortAlignments(int start,int end){
	int i,j,pivot;
	pivot=alignments_poses[start];
	i=start;j=end-1;
	while(i<j){
		while(i<j&&alignments[alignments_poses[j]].ref_pos>alignments[pivot].ref_pos)
			j--;
		if(i==j)
			break;
		alignments_poses[i]=alignments_poses[j];
		i++;
		while(i<j&&alignments[alignments_poses[i]].ref_pos<alignments[pivot].ref_pos)
			i++;
		if(i==j)
			break;
		alignments_poses[j]=alignments_poses[i];
		j--;
	}
	alignments_poses[i]=pivot;
	if(start<i)
		sortAlignments(start,i);
	if(i+1<end)
		sortAlignments(i+1,end);
}
int PEReads::lookupLowerBound(NaiveAlignmentFormat * alignments,int size,int ref_pos){
	int start,end,mid;
	start=-1;end=size;
	while(start<end-1){
		mid=(start+end)/2;
		if(ref_pos<=alignments[alignments_poses[mid]].ref_pos)
			end=mid;
		else
			start=mid;
	}
	return start;
}
int PEReads::lookupUpperBound(NaiveAlignmentFormat * alignments,int size,int ref_pos){
	int start,end,mid;
	start=-1;end=size;
	while(start<end-1){
		mid=(start+end)/2;
		if(ref_pos<alignments[alignments_poses[mid]].ref_pos)
			end=mid;
		else
			start=mid;
	}
	return end;
}
bool PEReads::findPositiveBestPair(){
	int i,start_pos,end_pos,read1_len,read2_len;
	int select_left,select_right;
	read1_len=Configure::read1.length();
	read2_len=Configure::read2.length();
	if(left_size<right_size){
		for(i=0;i!=left_size;i++)
			alignments_poses[i]=i;
		alignments=left_alignments;
		if(left_size>1)
			sortAlignments(0,left_size);
		for(i=0;i!=right_size;i++){
			start_pos=lookupLowerBound(left_alignments,left_size,right_alignments[i].ref_pos-high_bound+read2_len);
			end_pos=lookupUpperBound(left_alignments,left_size,right_alignments[i].ref_pos-low_bound+read2_len);
			if(start_pos+1>=end_pos)
				continue;
			select_left=start_pos+rand()%(end_pos-start_pos-1)+1;
			select_right=i;
			break;
		}
		if(i==right_size)
			return false;
		pe_alignment.start_poses[0]=left_alignments[select_left].ref_pos;
		pe_alignment.start_poses[1]=right_alignments[select_right].ref_pos;
		if(!pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len))
			return false;
		pe_alignment.insert_size=right_alignments[select_right].ref_pos-left_alignments[select_left].ref_pos+right_alignments[select_right].ref_read_len;
		pe_alignment.edits[0]=left_alignments[select_left].edit;
		pe_alignment.edits[1]=right_alignments[select_right].edit;
		pe_alignment.cigar[0]=left_alignments[select_left].cigar;
		pe_alignment.cigar[1]=right_alignments[select_right].cigar;
		return true;
	}
	for(i=0;i!=right_size;i++)
		alignments_poses[i]=i;
	alignments=right_alignments;
	if(right_size>1)
		sortAlignments(0,right_size);
	for(i=0;i!=left_size;i++){
		start_pos=lookupLowerBound(right_alignments,right_size,left_alignments[i].ref_pos+low_bound-read2_len);
		end_pos=lookupUpperBound(right_alignments,right_size,left_alignments[i].ref_pos+high_bound-read2_len);
		if(start_pos+1>=end_pos)
			continue;
		select_right=start_pos+rand()%(end_pos-start_pos-1)+1;
		select_left=i;
		break;
	}
	if(i==left_size)
		return false;
	pe_alignment.start_poses[0]=left_alignments[select_left].ref_pos;
	pe_alignment.start_poses[1]=right_alignments[select_right].ref_pos;
	if(!pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len))
		return false;
	pe_alignment.insert_size=right_alignments[select_right].ref_pos-left_alignments[select_left].ref_pos+right_alignments[select_right].ref_read_len;
	pe_alignment.edits[0]=left_alignments[select_left].edit;
	pe_alignment.edits[1]=right_alignments[select_right].edit;
	pe_alignment.cigar[0]=left_alignments[select_left].cigar;
	pe_alignment.cigar[1]=right_alignments[select_right].cigar;
	return true;
}
bool PEReads::findNegativeBestPair(){
	int i,start_pos,end_pos,read1_len,read2_len;
	int select_left,select_right;
	read1_len=Configure::read1.length();
	read2_len=Configure::read2.length();
	if(left_size<right_size){
		for(i=0;i!=left_size;i++)
			alignments_poses[i]=i;
		alignments=left_alignments;
		if(left_size>1)
			sortAlignments(0,left_size);
		for(i=0;i!=right_size;i++){
			start_pos=lookupLowerBound(left_alignments,left_size,right_alignments[i].ref_pos+low_bound-read1_len);
			end_pos=lookupUpperBound(left_alignments,left_size,right_alignments[i].ref_pos+high_bound-read1_len);
			if(start_pos+1>=end_pos)
				continue;
			select_left=start_pos+rand()%(end_pos-start_pos-1)+1;
			select_right=i;
			break;
		}
		if(i==right_size)
			return false;
		pe_alignment.start_poses[0]=left_alignments[select_left].ref_pos;
		pe_alignment.start_poses[1]=right_alignments[select_right].ref_pos;
		if(!pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len))
			return false;
		pe_alignment.insert_size=right_alignments[select_right].ref_pos-left_alignments[select_left].ref_pos-left_alignments[select_left].ref_read_len;
		pe_alignment.edits[0]=left_alignments[select_left].edit;
		pe_alignment.edits[1]=right_alignments[select_right].edit;
		pe_alignment.cigar[0]=left_alignments[select_left].cigar;
		pe_alignment.cigar[1]=right_alignments[select_right].cigar;
		return true;
	}
	for(i=0;i!=right_size;i++)
		alignments_poses[i]=i;
	alignments=right_alignments;
	if(right_size>1)
		sortAlignments(0,right_size);
	for(i=0;i!=left_size;i++){
		start_pos=lookupLowerBound(right_alignments,right_size,left_alignments[i].ref_pos-high_bound+read1_len);
		end_pos=lookupUpperBound(right_alignments,right_size,left_alignments[i].ref_pos-low_bound+read1_len);
		if(start_pos+1>=end_pos)
			continue;
		select_right=start_pos+rand()%(end_pos-start_pos-1)+1;
		select_left=i;
		break;
	}
	if(i==left_size)
		return false;
	pe_alignment.start_poses[0]=left_alignments[select_left].ref_pos;
	pe_alignment.start_poses[1]=right_alignments[select_right].ref_pos;
	if(!pRef_seq->countRefIndex(pe_alignment,read1_len,read2_len))
		return false;
	pe_alignment.insert_size=right_alignments[select_right].ref_pos-left_alignments[select_left].ref_pos-left_alignments[select_left].ref_read_len;
	pe_alignment.edits[0]=left_alignments[select_left].edit;
	pe_alignment.edits[1]=right_alignments[select_right].edit;
	pe_alignment.cigar[0]=left_alignments[select_left].cigar;
	pe_alignment.cigar[1]=right_alignments[select_right].cigar;
	return true;
}
void PEReads::printInserts(){
	int i;
	for(i=0;i!=inserts_size;i++){
		if(i&&i%10==0)
			cout<<endl;
		cout<<setw(8)<<left<<inserts[i]<<" ";
	}
	cout<<endl;
}
void PEReads::printNaiveAlignments(NaiveAlignmentFormat * alignments,int size){
	int i;
	for(i=0;i!=size;i++){
		cout<<alignments[i].ref_pos<<","<<alignments[i].edit<<","<<alignments[i].cigar<<endl;
	}
}
void PEReads::printPEAlignment(const PEAlignmentFormat & alignment){
	cout<<alignment.ref_indices[0]<<","<<alignment.start_poses[0]<<endl;
	cout<<alignment.cigar[0]<<endl;
	cout<<alignment.ref_indices[1]<<","<<alignment.start_poses[1]<<endl;
	cout<<alignment.cigar[1]<<endl;
}
