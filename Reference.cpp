/*
 * Reference::.cpp
 *
 *  Created on: Apr 21, 2014
 *      Author: laiwenqi
 */

#include "Reference.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
using namespace std;
Reference::Reference() {
	refs_start=new int[Configure::max_refs];
}

Reference::~Reference() {
}
char Reference::table[]={'A','T','G','C'};
char Reference::mask=0xDF;
bool Reference::purifyAndFillGap(){
	char c;
	int ref_size;
	string path,temp;
	ifstream input;
	path=Configure::dir+Configure::refs_file;
	input.open(path.c_str());
	if(!input) {
		cerr<<"Error:The reference file can't be open!"<<endl;
		return false;
	}
	srand(time(NULL));
	refs_start_size=0;
	total_N=ref_size=0;
	while(input>>c){
		switch(c){
		case 'A':
		case 'T':
		case 'G':
		case 'C':
			ref_str.push_back(c);
			ref_size++;
			break;
		case 'a':
		case 't':
		case 'g':
		case 'c':
			ref_str.push_back((char)(c&mask));
			ref_size++;
			break;
		case 'N':
			ref_str.push_back(table[rand()%4]);
			ref_size++;
			total_N++;
			break;
		case '>':
			refs_start[refs_start_size]=ref_size;
			refs_start_size++;
			getline(input,temp);
			break;
		default:
			break;
		}
	}
	input.close();
	refs_start[refs_start_size]=ref_size;
	refs_start_size++;
	return true;
}
int Reference::countRefIndexs(int * poses,const string & cigar,int * prefs,int len,int read_len){
	int i,j,k,delta,cigar_size,valid_size;
	valid_size=0;
	cigar_size=cigar.length();
	for(i=0;i!=len;i++){
		for(j=1;j!=refs_start_size;j++)
			if(refs_start[j]>poses[i])
				break;
		if(poses[i]+Configure::global_max_gap+read_len>refs_start[j]){
			delta=0;
			for(k=0;k!=cigar_size;k++){
				if(cigar[k]=='I')
					delta--;
				else if(cigar[k]=='D')
					delta++;
			}
			if(poses[i]+read_len+delta>refs_start[j])
				prefs[i]=-1;
			else
				prefs[i]=j-1;
		}
		else{
			prefs[i]=j-1;
		}
		if(prefs[i]!=-1){
			poses[i]-=refs_start[j-1];
			poses[i]++;
			valid_size++;
		}
		else{
			poses[i]=-1;
		}
	}
	return valid_size;
}
int Reference::countRefIndexs(int * poses,int * prefs,int len,int read_len){
	int i,j,valid_size;
	valid_size=0;
	for(i=0;i!=len;i++){
		for(j=1;j!=refs_start_size;j++)
			if(refs_start[j]>poses[i])
				break;
		if(poses[i]+read_len>refs_start[j])
			prefs[i]=-1;
		else
			prefs[i]=j-1;
		if(prefs[i]!=-1){
			poses[i]-=refs_start[j-1];
			poses[i]++;
			valid_size++;
		}
		else{
			poses[i]=-1;
		}
	}
	return valid_size;
}
int Reference::countRefIndexs(int * poses,int delta,int * refs,int len){
	int i,j,valid_size;
	valid_size=0;
	for(i=0;i!=len;i++){
		for(j=1;j!=refs_start_size;j++)
			if(refs_start[j]>poses[i])
				break;
		if(poses[i]+delta>refs_start[j])
			refs[i]=-1;
		else
			refs[i]=j-1;
		if(refs[i]!=-1){
			poses[i]-=refs_start[j-1];
			poses[i]++;
			valid_size++;
		}
		else{
			poses[i]=-1;
		}
	}
	return valid_size;
}
void Reference::countRefIndex(AlignmentFormat & alignment,int read_len){
	int i,delta,cigar_size,start_pos,ref_index;
	int start,mid,end;
	string cigar;
	start_pos=alignment.start_pos;
	start=0;end=refs_start_size;
	while(start<end-1){
		mid=(start+end)/2;
		if(start_pos<refs_start[mid])
			end=mid;
		else
			start=mid;
	}
	if(start_pos+Configure::global_max_gap+read_len>refs_start[end]){
		delta=0;
		cigar=alignment.cigar;
		cigar_size=cigar.length();
		for(i=0;i!=cigar_size;i++){
			if(cigar[i]=='I')
				delta--;
			else if(cigar[i]=='D')
				delta++;
		}
		if(start_pos+read_len+delta>refs_start[end])
			ref_index=-1;
		else
			ref_index=end-1;
	}
	else
		ref_index=end-1;
	if(ref_index!=-1){
		alignment.start_pos=start_pos-refs_start[ref_index];
		alignment.ref_index=ref_index;
	}
	else{
		alignment.ref_index=-1;
	}
}
void Reference::countRefIndex(int & ref_index,int & start_pos,const string & cigar,int read_len){
	int i,delta,cigar_size,pos,temp_ref_index;
	int start,mid,end;
	cigar_size=cigar.length();
	pos=start_pos;
	start=0;end=refs_start_size;
	while(start<end-1){
		mid=(start+end)/2;
		if(start_pos<refs_start[mid])
			end=mid;
		else
			start=mid;
	}
	if(pos+Configure::global_max_gap+read_len>refs_start[end]){
		delta=0;
		for(i=0;i!=cigar_size;i++){
			if(cigar[i]=='I')
				delta--;
			else if(cigar[i]=='D')
				delta++;
		}
		if(pos+read_len+delta>refs_start[end])
			temp_ref_index=-1;
		else
			temp_ref_index=end-1;
	}
	else
		temp_ref_index=end-1;
	if(temp_ref_index!=-1){
		start_pos=pos-refs_start[temp_ref_index];
		ref_index=temp_ref_index;
	}
	else{
		ref_index=-1;
	}
}
void Reference::countRefIndex(int & ref_index,int & start_pos,int read_len){
	int pos,temp_ref_index;
	int start,end,mid;
	pos=start_pos;
	start=0;end=refs_start_size;
	while(start<end-1){
		mid=(start+end)/2;
		if(start_pos<refs_start[mid])
			end=mid;
		else
			start=mid;
	}
	if(pos+Configure::global_max_gap+read_len>refs_start[end])
		temp_ref_index=-1;
	else
		temp_ref_index=end-1;
	if(temp_ref_index!=-1){
		start_pos=pos-refs_start[temp_ref_index];
		ref_index=temp_ref_index;
	}
	else{
		ref_index=-1;
	}
}
void Reference::printRefsStart(){
	int i;
	for(i=0;i!=refs_start_size;i++)
		cout<<refs_start[i]<<" ";
	cout<<endl;
}
bool Reference::countRefIndex(PEAlignmentFormat & pe_alignment,int read1_len, int read2_len){
	int i,delta,cigar_size,start_pos,ref_index;
	int start,end,mid;
	string cigar;
	start_pos=pe_alignment.start_poses[0];
	start=0;end=refs_start_size;
	while(start<end-1){
		mid=(start+end)/2;
		if(start_pos<refs_start[mid])
			end=mid;
		else
			start=mid;
	}
	if(start_pos+Configure::global_max_gap+read1_len>refs_start[end]){
		delta=0;
		cigar=pe_alignment.cigar[0];
		cigar_size=cigar.length();
		for(i=0;i!=cigar_size;i++){
			if(cigar[i]=='I')
				delta--;
			else if(cigar[i]=='D')
				delta++;
		}
		if(start_pos+read1_len+delta>refs_start[end])
			ref_index=-1;
		else
			ref_index=end-1;
	}
	else
		ref_index=end-1;
	if(ref_index!=-1){
		pe_alignment.start_poses[0]=start_pos-refs_start[ref_index];
		pe_alignment.ref_indices[0]=ref_index;
	}
	else{
		pe_alignment.ref_indices[0]=-1;
	}
	if(ref_index==-1)
		return false;
	start_pos=pe_alignment.start_poses[1];
	start=0;end=refs_start_size;
	while(start<end-1){
		mid=(start+end)/2;
		if(start_pos<refs_start[mid])
			end=mid;
		else
			start=mid;
	}
	if(start_pos+Configure::global_max_gap+read2_len>refs_start[end]){
		delta=0;
		cigar=pe_alignment.cigar[1];
		cigar_size=cigar.length();
		for(i=0;i!=cigar_size;i++){
			if(cigar[i]=='I')
				delta--;
			else if(cigar[i]=='D')
				delta++;
		}
		if(start_pos+read2_len+delta>refs_start[end])
			ref_index=-1;
		else
			ref_index=end-1;
	}
	else
		ref_index=end-1;
	if(ref_index!=-1){
		pe_alignment.start_poses[1]=start_pos-refs_start[ref_index];
		pe_alignment.ref_indices[1]=ref_index;
	}
	else{
		pe_alignment.ref_indices[1]=-1;
	}
	if(ref_index!=pe_alignment.ref_indices[0])
		return false;
	return true;
}
