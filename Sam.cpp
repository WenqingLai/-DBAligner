/*
 * Sam.cpp
 *
 *  Created on: Mar 19, 2014
 *      Author: laiwenqi
 */

#include "Sam.h"
#include "Configure.h"
#include<sstream>
#include<iostream>
int Sam::unmapped=0x4;
int Sam::mapped_neg=0x10;
int Sam::pos_secondary=0x100;
int Sam::neg_secondary=0x110;
int Sam::pe_positive_mapped_f=0x63;
int Sam::pe_positive_mapped_s=0x93;
int Sam::pe_negative_mapped_f=0x53;
int Sam::pe_negative_mapped_s=0xA3;
int Sam::pe_disc_ff_f=0x41;
int Sam::pe_disc_fr_f=0x61;
int Sam::pe_disc_rf_f=0x51;
int Sam::pe_disc_rr_f=0x71;
int Sam::pe_disc_ff_s=0x81;
int Sam::pe_disc_fr_s=0x91;
int Sam::pe_disc_rf_s=0xA1;
int Sam::pe_disc_rr_s=0xB1;

Sam::Sam() {
	// TODO Auto-generated constructor stub
}

Sam::~Sam() {
	// TODO Auto-generated destructor stub
}
bool Sam::writeOpen(const char * fn){
	output.open(fn);
	if(!output)
		return false;
	return true;
}
void Sam::writeHeader(int ref_lens[],int size){
	int i;
	for(i=1;i!=size;i++)
		output<<"@SQ"<<"\tSN:"<<(i-1)<<"\tLN:"<<ref_lens[i]<<endl;
}
void Sam::writeUnmappedAlignment(int read_num,const string & read){
	output<<"r"<<read_num<<"\t"<<unmapped;
	output<<"\t*"<<"\t0"<<"\t0\t*"<<"\t*\t0";
	output<<"\t0\t"<<read<<"\t*"<<endl;
}
void Sam::writePosMappedAlignment(int read_num,int edit,const string & cigar,const string & read,int * prefs,int * poses,int size){
	int i;
	output<<"r"<<read_num<<"\t0";
	for(i=0;i!=size;i++){
		if(prefs[i]!=-1)
			break;
	}
	output<<"\t"<<prefs[i]<<"\t"<<poses[i];
	output<<"\t0\t"<<compressCigar(cigar)<<"\t*\t0";
	output<<"\t0\t"<<read<<"\t*";
	output<<"\tNM:i:"<<edit;
	output<<"\tXI:B:s";
	for(i=0;i!=size;i++){
		if(prefs[i]==-1)
			continue;
		output<<","<<prefs[i];
	}
	output<<"\tXS:B:i";
	for(i=0;i!=size;i++){
		if(poses[i]==-1)
			continue;
		output<<","<<poses[i];
	}
	output<<endl;
}
void Sam::writeNegMappedAlignment(int read_num,int edit,const string & cigar,const string & read,int * prefs,int * poses,int size){
	int i;
	output<<"r"<<read_num<<"\t"<<mapped_neg;
	for(i=0;i!=size;i++){
		if(prefs[i]!=-1)
			break;
	}
	output<<"\t"<<prefs[i]<<"\t"<<poses[i];
	output<<"\t0\t"<<compressCigar(cigar)<<"\t*\t0";
	output<<"\t0\t"<<read<<"\t*";
	output<<"\tNM:i:"<<edit;
	output<<"\tXI:B:s";
	for(i=0;i!=size;i++){
		if(prefs[i]==-1)
			continue;
		output<<","<<prefs[i];
	}
	output<<"\tXS:B:i";
	for(i=0;i!=size;i++){
		if(poses[i]==-1)
			continue;
		output<<","<<poses[i];
	}
	output<<endl;
}
void Sam::writePosPerfectAlignment(int read_num,const string & read,int * prefs,int * poses,int size){
	int i;
	output<<"r"<<read_num<<"\t0";
	for(i=0;i!=size;i++){
		if(prefs[i]!=-1)
			break;
	}
	output<<"\t"<<prefs[i]<<"\t"<<poses[i];
	output<<"\t0\t"<<read.length()<<"=\t*\t0";
	output<<"\t0\t"<<read<<"\t*";
	output<<"\tNM:i:0";
	output<<"\tXI:B:s";
	for(i=0;i!=size;i++){
		if(prefs[i]==-1)
			continue;
		output<<","<<prefs[i];
	}
	output<<"\tXS:B:i";
	for(i=0;i!=size;i++){
		if(poses[i]==-1)
			continue;
		output<<","<<poses[i];
	}
	output<<endl;
}
void Sam::writeNegPerfectAlignment(int read_num,const string & read,int * prefs,int * poses,int size){
	int i;
	output<<"r"<<read_num<<"\t"<<mapped_neg;
	for(i=0;i!=size;i++){
		if(prefs[i]!=-1)
			break;
	}
	output<<"\t"<<prefs[i]<<"\t"<<poses[i];
	output<<"\t0\t"<<read.length()<<"=\t*\t0";
	output<<"\t0\t"<<read<<"\t*";
	output<<"\tNM:i:0";
	output<<"\tXI:B:s";
	for(i=0;i!=size;i++){
		if(prefs[i]==-1)
			continue;
		output<<","<<prefs[i];
	}
	output<<"\tXS:B:i";
	for(i=0;i!=size;i++){
		if(poses[i]==-1)
			continue;
		output<<","<<poses[i];
	}
	output<<endl;
}
void Sam::writePosMappedAlignment(int read_num,const string & read,const AlignmentFormat & alignment){
	output<<"r"<<read_num<<"\t0";
	output<<"\t"<<alignment.ref_index<<"\t"<<alignment.start_pos+1;
	output<<"\t0\t"<<compressCigar(alignment.cigar)<<"\t*\t0";
	output<<"\t0\t"<<read<<"\t*";
	output<<"\tNM:i:"<<alignment.edit<<endl;
}
void Sam::writeNegMappedAlignment(int read_num,const string & read,const AlignmentFormat & alignment){
	output<<"r"<<read_num<<"\t"<<mapped_neg;
	output<<"\t"<<alignment.ref_index<<"\t"<<alignment.start_pos+1;
	output<<"\t0\t"<<compressCigar(alignment.cigar)<<"\t*\t0";
	output<<"\t0\t"<<read<<"\t*";
	output<<"\tNM:i:"<<alignment.edit<<endl;
}
void Sam::writePosSecMappedAlignment(int read_num,const AlignmentFormat & alignment){
	output<<"r"<<read_num<<"\t"<<pos_secondary;
	output<<"\t"<<alignment.ref_index<<"\t"<<alignment.start_pos+1;
	output<<"\t0\t"<<compressCigar(alignment.cigar)<<"\t*\t0";
	output<<"\t0\t*\t*\tNM:i:"<<alignment.edit<<endl;
}
void Sam::writeNegSecMappedAlignment(int read_num,const AlignmentFormat & alignment){
	output<<"r"<<read_num<<"\t"<<neg_secondary;
	output<<"\t"<<alignment.ref_index<<"\t"<<alignment.start_pos+1;
	output<<"\t0\t"<<compressCigar(alignment.cigar)<<"\t*\t0";
	output<<"\t0\t*\t*\tNM:i:"<<alignment.edit<<endl;
}
void Sam::writePosPEAlignment(int read_num,const PEAlignmentFormat & pe_alignment){
	output<<"r"<<read_num<<"_1\t";
	output<<pe_positive_mapped_f;
	output<<"\t"<<pe_alignment.ref_indices[0]<<"\t"<<pe_alignment.start_poses[0]+1;
	output<<"\t0\t"<<compressCigar(pe_alignment.cigar[0]);
	output<<"\t=\t"<<pe_alignment.start_poses[1]+1;
	output<<"\t"<<pe_alignment.insert_size<<"\t"<<Configure::read1<<"\t*";
	output<<"\tNM:i:"<<pe_alignment.edits[0]<<endl;
	output<<"r"<<read_num<<"_2\t";
	output<<pe_positive_mapped_s;
	output<<"\t"<<pe_alignment.ref_indices[1]<<"\t"<<pe_alignment.start_poses[1]+1;
	output<<"\t0\t"<<compressCigar(pe_alignment.cigar[1]);
	output<<"\t=\t"<<pe_alignment.start_poses[0]+1;
	output<<"\t"<<-pe_alignment.insert_size<<"\t"<<Configure::rev_read2<<"\t*";
	output<<"\tNM:i:"<<pe_alignment.edits[1]<<endl;
}
void Sam::writeNegPEAlignment(int read_num,const PEAlignmentFormat & pe_alignment){
	output<<"r"<<read_num<<"_1\t";
	output<<pe_negative_mapped_f;
	output<<"\t"<<pe_alignment.ref_indices[0]<<"\t"<<pe_alignment.start_poses[0]+1;
	output<<"\t0\t"<<compressCigar(pe_alignment.cigar[0]);
	output<<"\t=\t"<<pe_alignment.start_poses[1]+1;
	output<<"\t"<<pe_alignment.insert_size<<"\t"<<Configure::rev_read1<<"\t*";
	output<<"\tNM:i:"<<pe_alignment.edits[0]<<endl;
	output<<"r"<<read_num<<"_2\t";
	output<<pe_negative_mapped_s;
	output<<"\t"<<pe_alignment.ref_indices[1]<<"\t"<<pe_alignment.start_poses[1]+1;
	output<<"\t0\t"<<compressCigar(pe_alignment.cigar[1]);
	output<<"\t=\t"<<pe_alignment.start_poses[0]+1;
	output<<"\t"<<-pe_alignment.insert_size<<"\t"<<Configure::read2<<"\t*";
	output<<"\tNM:i:"<<pe_alignment.edits[1]<<endl;
}
void Sam::writePEDiscordance(int read_num,const PEAlignmentFormat & pe_alignment,bool r1_positive,bool r2_positive){
	output<<"r"<<read_num<<"_1\t";
	if(r1_positive){
		if(r2_positive)
			output<<pe_disc_ff_f;
		else
			output<<pe_disc_fr_f;
	}
	else{
		if(r2_positive)
			output<<pe_disc_rf_f;
		else
			output<<pe_disc_rr_f;
	}
	output<<"\t"<<pe_alignment.ref_indices[0]<<"\t"<<pe_alignment.start_poses[0]+1;
	output<<"\t0\t"<<compressCigar(pe_alignment.cigar[0]);
	output<<"\t"<<pe_alignment.ref_indices[1]<<"\t"<<pe_alignment.start_poses[1]+1;
	output<<"\t0\t";
	if(r1_positive)
		output<<Configure::read1<<"\t*";
	else
		output<<Configure::rev_read1<<"\t*";
	output<<"\tNM:i:"<<pe_alignment.edits[0]<<endl;
	output<<"r"<<read_num<<"_2\t";
	if(r1_positive){
		if(r2_positive)
			output<<pe_disc_ff_s;
		else
			output<<pe_disc_fr_s;
	}
	else{
		if(r2_positive)
			output<<pe_disc_rf_s;
		else
			output<<pe_disc_rr_s;
	}
	output<<"\t"<<pe_alignment.ref_indices[1]<<"\t"<<pe_alignment.start_poses[1]+1;
	output<<"\t0\t"<<compressCigar(pe_alignment.cigar[1]);
	output<<"\t"<<pe_alignment.ref_indices[0]<<"\t"<<pe_alignment.start_poses[0]+1;
	output<<"\t0\t";
	if(r2_positive)
		output<<Configure::read2<<"\t*";
	else
		output<<Configure::rev_read2<<"\t*";
	output<<"\tNM:i:"<<pe_alignment.edits[1]<<endl;
}
void Sam::writePEUnmappedAlignment(int read_num,const PEAlignmentFormat & pe_alignment,bool r1_positive,bool r2_positive){
	if(pe_alignment.ref_indices[0]==-1){
		if(pe_alignment.ref_indices[1]==-1){
			output<<"r"<<read_num<<"_1\t";
			output<<0x4D<<"\t*\t0\t0\t*\t*\t0\t0\t"<<Configure::read1<<"\t*"<<endl;
			output<<"r"<<read_num<<"_2\t";
			output<<0x8D<<"\t*\t0\t0\t*\t*\t0\t0\t"<<Configure::read2<<"\t*"<<endl;
		}
		else if(r2_positive){
			output<<"r"<<read_num<<"_1\t";
			output<<0x45<<"\t*\t0\t0\t*\t"<<pe_alignment.ref_indices[1]<<"\t";
			output<<pe_alignment.start_poses[1]+1<<"\t0\t"<<Configure::read1<<"\t*"<<endl;
			output<<"r"<<read_num<<"_2\t";
			output<<0x89<<"\t"<<pe_alignment.ref_indices[1]<<"\t"<<pe_alignment.start_poses[1]+1;
			output<<"\t0\t"<<compressCigar(pe_alignment.cigar[1]);
			output<<"\t*\t0\t0\t"<<Configure::read2<<"\t*"<<"\tNM:i:"<<pe_alignment.edits[1];
			output<<endl;
		}
		else{
			output<<"r"<<read_num<<"_1\t";
			output<<0x65<<"\t*\t0\t0\t*\t"<<pe_alignment.ref_indices[1]<<"\t";
			output<<pe_alignment.start_poses[1]+1<<"\t0\t"<<Configure::read1<<"\t*"<<endl;
			output<<"r"<<read_num<<"_2\t";
			output<<0x99<<"\t"<<pe_alignment.ref_indices[1]<<"\t"<<pe_alignment.start_poses[1]+1;
			output<<"\t0\t"<<compressCigar(pe_alignment.cigar[1]);
			output<<"\t*\t0\t0\t"<<Configure::rev_read2<<"\t*"<<"\tNM:i:"<<pe_alignment.edits[1];
			output<<endl;
		}
	}
	else if(pe_alignment.ref_indices[1]==-1){
		if(r1_positive){
			output<<"r"<<read_num<<"_1\t";
			output<<0x49<<"\t"<<pe_alignment.ref_indices[0]<<"\t"<<pe_alignment.start_poses[0]+1;
			output<<"\t0\t"<<compressCigar(pe_alignment.cigar[0]);
			output<<"\t*\t0\t0\t"<<Configure::read1<<"\t*";
			output<<"\tNM:i:"<<pe_alignment.edits[0]<<endl;
			output<<"r"<<read_num<<"_2\t";
			output<<0x85<<"\t*\t0\t0\t*\t"<<pe_alignment.ref_indices[0]<<"\t";
			output<<pe_alignment.start_poses[0]+1<<"\t0\t"<<Configure::read2<<"\t*"<<endl;
		}
		else{
			output<<"r"<<read_num<<"_1\t";
			output<<0x59<<"\t"<<pe_alignment.ref_indices[0]<<"\t"<<pe_alignment.start_poses[0]+1;
			output<<"\t0\t"<<compressCigar(pe_alignment.cigar[0]);
			output<<"\t*\t0\t0\t"<<Configure::rev_read1<<"\t*";
			output<<"\tNM:i:"<<pe_alignment.edits[0]<<endl;
			output<<"r"<<read_num<<"_2\t";
			output<<0xA5<<"\t*\t0\t0\t*\t"<<pe_alignment.ref_indices[0]<<"\t";
			output<<pe_alignment.start_poses[0]+1<<"\t0\t"<<Configure::read2<<"\t*"<<endl;
		}
	}
}
void Sam::writePEUnmappedAlignment(int read_num){
	output<<"r"<<read_num<<"_1\t";
	output<<0x4D<<"\t*\t0\t0\t*\t*\t0\t0\t"<<Configure::read1<<"\t*"<<endl;
	output<<"r"<<read_num<<"_2\t";
	output<<0x8D<<"\t*\t0\t0\t*\t*\t0\t0\t"<<Configure::read2<<"\t*"<<endl;
}
void Sam::writeClose(){
	output.close();
}
string Sam::compressCigar(const string & raw_cigar){
	char prev;
	int i,count,len;
	stringstream output;
	const char * praw_cigar;
	praw_cigar=raw_cigar.c_str();
	len=raw_cigar.length();
	prev=praw_cigar[0];
	count=1;
	for(i=1;i!=len;i++){
		if(praw_cigar[i]==prev){
			count++;
		}
		else{
			output<<count<<prev;
			count=1;
		}
		prev=praw_cigar[i];
	}
	output<<count<<prev;
	return output.str();
}
string Sam::compressCigar(const string * raw_cigar){
	char prev;
	int i,count,len;
	stringstream output;
	const char * praw_cigar;
	praw_cigar=raw_cigar->c_str();
	len=raw_cigar->length();
	prev=praw_cigar[0];
	count=1;
	for(i=1;i!=len;i++){
		if(praw_cigar[i]==prev){
			count++;
		}
		else{
			output<<count<<prev;
			count=1;
		}
		prev=praw_cigar[i];
	}
	output<<count<<prev;
	return output.str();
}

