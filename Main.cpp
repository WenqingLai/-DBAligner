/*
 * Main.cpp
 *
 *  Created on: Apr 21, 2014
 *      Author: laiwenqi
 */
#include <iostream>
#include <fstream>
#include <sys/resource.h>
#include <sys/time.h>
#include "Reference.h"
#include "DNA.h"
#include "Sam.h"
#include "Configure.h"
#include "SubRead.h"
#include "Graph.h"
#include "HashTable.h"
#include "PEReads.h"
using namespace std;
double cputime();
double realtime();
void alignSingleReads();
void alignSingleReads_All();
void alignPairedReads();
void printMarks(bool marks[],int len);
string blank_line;
int main(int argc,char * argv[]){
	ifstream input;
	string read,rev_read,temp,path;
	double t_real;
	t_real=realtime();
	blank_line.append(70,' ');
	if(argc<=1){
		if(!Configure::initialize())
			return 0;
	}
	else{
		if(!Configure::initialize(argv[1]))
			return 0;
	}
	if(Configure::is_pe){
		Configure::printPEConfiguration();
		alignPairedReads();
	}
	else if(Configure::report_all){
		Configure::printSingleConfiguration();
		alignSingleReads_All();
	}
	else{
		Configure::printSingleConfiguration();
		alignSingleReads();
	}
	cout<<"real time:"<<realtime()-t_real<<",CPU:"<<cputime()<<endl;
	return 0;
}
void alignSingleReads(){
	string path,temp,read,rev_read;
	ifstream input;
	int total,read_len,mapped_num;
	int i,count1,count2,temp_edit,temp_gap;
	bool * temp_marks1,* temp_marks2,mapped;
	int temp_marks_len;
	AlignmentFormat temp_alignment;
	Reference * pRef_seq;
	Sam output;
	SubRead * pSub_read;
	Graph * pGraph;
	HashTable * pHash_table;
	path=Configure::dir+Configure::refs_file;
	input.open(path.c_str());
	if(!input){
		cout<<"The reference file can't be opened!"<<endl;
		return;
	}
	else
		input.close();
	path=Configure::dir+Configure::reads_file;
	input.open(path.c_str());
	if(!input){
		cout<<"can't open reads file."<<endl;
		return;
	}
	pRef_seq=new Reference();
	path=Configure::dir+Configure::align_file;
	if(!output.writeOpen(path.c_str())){
		cout<<"The output file can't be opened!"<<endl;
		return;
	}
	if(!pRef_seq->purifyAndFillGap())
		return;
	pHash_table=new HashTable(pRef_seq);
	pGraph=new Graph(pRef_seq,pHash_table);
	pSub_read=new SubRead(pRef_seq,pHash_table,pGraph);
	temp_marks1=new bool[Configure::max_read_len];
	temp_marks2=new bool[Configure::max_read_len];
	total=mapped_num=0;
	output.writeHeader(pRef_seq->refs_start,pRef_seq->refs_start_size);
	while(getline(input,temp)){
		getline(input,read);
		if(Configure::reads_fmt==FileFormat::fq){
			getline(input,temp);
			getline(input,temp);
		}
		total++;
		read_len=read.length();
		if(read_len<Configure::seed_size){
			output.writeUnmappedAlignment(total,read);
			continue;
		}
		Configure::global_max_edit=Configure::global_edit_ratio*read_len;
		Configure::global_max_gap=Configure::global_gap_ratio*read_len;
		if(total%10000==0){
			cout<<" total:"<<total<<",mapped:"<<mapped_num<<"\r";
			cout.flush();
		}
		rev_read=DNA::reverseComplete(read);
		temp_marks_len=read_len-Configure::seed_size+1;
		pGraph->makeMarks(read,temp_marks1,temp_marks_len);
		pGraph->makeMarks(rev_read,temp_marks2,temp_marks_len);
		count1=count2=0;
		for(i=0;i!=temp_marks_len;i++){
			if(temp_marks1[i])
				count1++;
			if(temp_marks2[i])
				count2++;
		}
		mapped=false;
		if(count1>count2){
			if(pGraph->searchGraph(read,temp_marks1,temp_marks_len)&&pGraph->resolveAlignment()){
				mapped=true;
				temp_alignment=pGraph->alignment;
			}
			temp_edit=Configure::global_max_edit;
			temp_gap=Configure::global_max_gap;
			if(mapped){
				Configure::global_max_edit=temp_alignment.edit;
				if(Configure::global_max_gap>temp_alignment.edit){
					Configure::global_max_gap=temp_alignment.edit;
				}
			}
			if(pGraph->searchGraph(rev_read,temp_marks2,temp_marks_len)&&pGraph->resolveAlignment()){
				mapped_num++;
				output.writeNegMappedAlignment(total,rev_read,pGraph->alignment);
				continue;
			}
			if(mapped){
				mapped_num++;
				output.writePosMappedAlignment(total,read,temp_alignment);
				continue;
			}
			Configure::global_max_edit=temp_edit;
			Configure::global_max_gap=temp_gap;
			if(pSub_read->alignRead(read,temp_marks1,temp_marks_len)){
				mapped_num++;
				output.writePosMappedAlignment(total,read,pSub_read->alignment);
				continue;
			}
			if(pSub_read->alignRead(rev_read,temp_marks2,temp_marks_len)){
				mapped_num++;
				output.writeNegMappedAlignment(total,rev_read,pSub_read->alignment);
				continue;
			}
		}
		else{
			if(pGraph->searchGraph(rev_read,temp_marks2,temp_marks_len)&&pGraph->resolveAlignment()){
				mapped=true;
				temp_alignment=pGraph->alignment;
			}
			temp_edit=Configure::global_max_edit;
			temp_gap=Configure::global_max_gap;
			if(mapped){
				Configure::global_max_edit=temp_alignment.edit;
				if(Configure::global_max_gap>temp_alignment.edit)
					Configure::global_max_gap=temp_alignment.edit;
			}
			if(pGraph->searchGraph(read,temp_marks1,temp_marks_len)&&pGraph->resolveAlignment()){
				mapped_num++;
				output.writePosMappedAlignment(total,read,pGraph->alignment);
				continue;
			}
			if(mapped){
				mapped_num++;
				output.writeNegMappedAlignment(total,rev_read,temp_alignment);
				continue;
			}
			Configure::global_max_edit=temp_edit;
			Configure::global_max_gap=temp_gap;
			if(pSub_read->alignRead(read,temp_marks1,temp_marks_len)){
				mapped_num++;
				output.writePosMappedAlignment(total,read,pSub_read->alignment);
				continue;
			}
			if(pSub_read->alignRead(rev_read,temp_marks2,temp_marks_len)){
				mapped_num++;
				output.writeNegMappedAlignment(total,rev_read,pSub_read->alignment);
				continue;
			}
		}
		output.writeUnmappedAlignment(total,read);
	}
	input.close();
	output.writeClose();
	cout.precision(2);
	cout<<fixed;
	cout<<blank_line<<endl;
	cout<<"total:"<<total<<",mapped:"<<mapped_num<<",mapped_ratio:"<<((float)mapped_num/total)*100<<"%"<<endl;
}
void alignPairedReads(){
	string path,temp,read1,read2;
	ifstream input1,input2;
	int total,mapped_num,disc_num,chimeric_num;
	AlignedType type;
	Reference * pRef_seq;
	Sam output;
	SubRead * pSub_read;
	Graph * pGraph;
	HashTable * pHash_table;
	PEReads * pPe_reads;
	bool estimate_not_enough;
	path=Configure::dir+Configure::refs_file;
	input1.open(path.c_str());
	if(!input1){
		cout<<"The reference file can't be opened!"<<endl;
		return;
	}
	else
		input1.close();
	path=Configure::dir+Configure::left_reads;
	input1.open(path.c_str());
	if(!input1){
		cout<<"can't open reads file 1."<<endl;
		return;
	}
	path=Configure::dir+Configure::right_reads;
	input2.open(path.c_str());
	if(!input2){
		cout<<"can't open reads file 2."<<endl;
		return;
	}
	pRef_seq=new Reference();
	path=Configure::dir+Configure::align_file;
	if(!output.writeOpen(path.c_str())){
		cout<<"The output file can't be opened!"<<endl;
		return;
	}
	if(!pRef_seq->purifyAndFillGap())
		return;
	pHash_table=new HashTable(pRef_seq);
	pGraph=new Graph(pRef_seq,pHash_table);
	pSub_read=new SubRead(pRef_seq,pHash_table, pGraph);
	pPe_reads=new PEReads(pRef_seq,pGraph,pSub_read);
	total=mapped_num=0;
	estimate_not_enough=true;
	if(!Configure::bound_designated){
		while(getline(input1,temp)&&getline(input2,temp)){
			getline(input1,Configure::read1);getline(input2,Configure::read2);
			if(Configure::reads_fmt==FileFormat::fq){
				getline(input1,temp);getline(input2,temp);
				getline(input1,temp);getline(input2,temp);
			}
			if(Configure::read1.length()<(unsigned int)Configure::seed_size||Configure::read2.length()<(unsigned int)Configure::seed_size)
				continue;
			Configure::rev_read1=DNA::reverseComplete(Configure::read1);
			Configure::rev_read2=DNA::reverseComplete(Configure::read2);
			if(pPe_reads->estimateInsert()){
				estimate_not_enough=false;
				break;
			}
		}
		if(estimate_not_enough)
			pPe_reads->caculateInsert();
		cout<<"average insert:"<<pPe_reads->average<<",standard deviation:"<<pPe_reads->st_deviation<<endl;
	}
	else{
		pPe_reads->low_bound=Configure::low_bound;
		pPe_reads->high_bound=Configure::high_bound;
	}
	input1.clear();
	input2.clear();
	input1.seekg(0,ios::beg);
	input2.seekg(0,ios::beg);
	mapped_num=chimeric_num=disc_num=total=0;
	output.writeHeader(pRef_seq->refs_start,pRef_seq->refs_start_size);
	while(getline(input1,temp)&&getline(input2,temp)){
		getline(input1,Configure::read1);getline(input2,Configure::read2);
		if(Configure::reads_fmt==FileFormat::fq){
			getline(input1,temp);getline(input2,temp);
			getline(input1,temp);getline(input2,temp);
		}
		if(Configure::read1.length()<(unsigned int)Configure::seed_size||Configure::read2.length()<(unsigned int)Configure::seed_size){
			output.writePEUnmappedAlignment(total);
			continue;
		}
		Configure::rev_read1=DNA::reverseComplete(Configure::read1);
		Configure::rev_read2=DNA::reverseComplete(Configure::read2);
		total++;
		if(total%10000==0){
			cout<<" total:"<<total<<",mapped:"<<mapped_num<<",discordant:";
			cout<<disc_num<<",chimeric:"<<chimeric_num<<"\r";
			cout.flush();
		}
		type=pPe_reads->alignPERead();
		switch(type){
		case pe_mapped:
			mapped_num++;
			if(pPe_reads->pe_alignment.insert_size>0)
				output.writePosPEAlignment(total,pPe_reads->pe_alignment);
			else
				output.writeNegPEAlignment(total,pPe_reads->pe_alignment);
			break;
		case pe_discordant:
			disc_num++;
			output.writePEDiscordance(total,pPe_reads->pe_alignment,pPe_reads->left_positive,pPe_reads->right_positive);
			break;
		case pe_chimeric:
			chimeric_num++;
			if(pPe_reads->pe_alignment.insert_size>0)
				output.writePosPEAlignment(total,pPe_reads->pe_alignment);
			else
				output.writeNegPEAlignment(total,pPe_reads->pe_alignment);
			break;
		default:
			output.writePEUnmappedAlignment(total,pPe_reads->pe_alignment,pPe_reads->left_positive,pPe_reads->right_positive);
			break;
		};
	}
	input1.close();
	input2.close();
	output.writeClose();
	cout.precision(2);
	cout<<fixed;
	cout<<blank_line<<endl;
	cout<<"total:"<<total<<",mapped:"<<mapped_num<<",discordant:"<<disc_num<<",chimeric:"<<chimeric_num<<endl;
	cout<<"mapped_ratio:"<<((float)mapped_num/total)*100<<"%";
	cout<<",discordant_ratio:"<<((float)disc_num/total)*100<<"%";
	cout<<",chimeric_ratio:"<<((float)chimeric_num/total)*100<<"%"<<endl;
}
void alignSingleReads_All(){
	string path,temp,read,rev_read;
	ifstream input;
	int total,read_len,mapped_num,aligned_poses;
	int i,alignments_size,temp_marks_len;
	bool * pmarks,* nmarks,mapped;
	Reference * pRef_seq;
	Sam output;
	SubRead * pSub_read;
	Graph * pGraph;
	HashTable * pHash_table;
	NaiveAlignmentFormat * alignments;
	AlignmentFormat alignment;
	alignments=new NaiveAlignmentFormat[Configure::max_local_results];
	path=Configure::dir+Configure::refs_file;
	input.open(path.c_str());
	if(!input){
		cout<<"The reference file can't be opened!"<<endl;
		return;
	}
	else
		input.close();
	path=Configure::dir+Configure::reads_file;
	input.open(path.c_str());
	if(!input){
		cout<<"can't open reads file."<<endl;
		return;
	}
	pRef_seq=new Reference();
	path=Configure::dir+Configure::align_file;
	if(!output.writeOpen(path.c_str())){
		cout<<"The output file can't be opened!"<<endl;
		return;
	}
	if(!pRef_seq->purifyAndFillGap())
		return;
	pHash_table=new HashTable(pRef_seq);
	pGraph=new Graph(pRef_seq,pHash_table);
	pSub_read=new SubRead(pRef_seq,pHash_table,pGraph);
	pmarks=new bool[Configure::max_read_len];
	nmarks=new bool[Configure::max_read_len];
	total=mapped_num=aligned_poses=0;
	output.writeHeader(pRef_seq->refs_start,pRef_seq->refs_start_size);
	while(getline(input,temp)){
		getline(input,read);
		if(Configure::reads_fmt==FileFormat::fq){
			getline(input,temp);
			getline(input,temp);
		}
		total++;
		if(read_len<Configure::seed_size){
			output.writeUnmappedAlignment(total,read);
			continue;
		}
		read_len=read.length();
		Configure::global_max_edit=Configure::global_edit_ratio*read_len;
		Configure::global_max_gap=Configure::global_gap_ratio*read_len;
		if(total%10000==0){
			cout<<" total reads:"<<total<<",mapped reads:"<<mapped_num<<",mapped positions:"<<aligned_poses<<"\r";
			cout.flush();
		}
		mapped=false;
		temp_marks_len=read_len-Configure::seed_size+1;
		rev_read=DNA::reverseComplete(read);
		pGraph->makeMarks(read,pmarks,temp_marks_len);
		pGraph->makeMarks(rev_read,nmarks,temp_marks_len);
		if(pGraph->searchGraph(read,pmarks,temp_marks_len)){
			alignments_size=pGraph->resolveAlignments(alignments);
			for(i=0;i!=alignments_size;i++){
				alignment.start_pos=alignments[i].ref_pos;
				alignment.cigar=alignments[i].cigar;
				alignment.edit=alignments[i].edit;
				pRef_seq->countRefIndex(alignment,read_len);
				if(alignment.ref_index!=-1){
					output.writePosMappedAlignment(total,read,alignment);
					mapped=true;
					aligned_poses++;
					i++;
					break;
				}
			}
			for(;i!=alignments_size;i++){
				alignment.start_pos=alignments[i].ref_pos;
				alignment.cigar=alignments[i].cigar;
				alignment.edit=alignments[i].edit;
				pRef_seq->countRefIndex(alignment,read_len);
				if(alignment.ref_index!=-1){
					output.writePosSecMappedAlignment(total,alignment);
					aligned_poses++;
				}
			}
		}
		if(pGraph->searchGraph(rev_read,nmarks,temp_marks_len)){
			alignments_size=pGraph->resolveAlignments(alignments);
			for(i=0;i!=alignments_size;i++){
				alignment.start_pos=alignments[i].ref_pos;
				alignment.cigar=alignments[i].cigar;
				alignment.edit=alignments[i].edit;
				pRef_seq->countRefIndex(alignment,read_len);
				if(alignment.ref_index!=-1){
					if(!mapped){
						output.writeNegMappedAlignment(total,rev_read,alignment);
						mapped=true;
						aligned_poses++;
						i++;
						break;
					}
					else{
						output.writeNegSecMappedAlignment(total,alignment);
						mapped=true;
						aligned_poses++;
						i++;
						break;
					}
				}
			}
			for(;i!=alignments_size;i++){
				alignment.start_pos=alignments[i].ref_pos;
				alignment.cigar=alignments[i].cigar;
				alignment.edit=alignments[i].edit;
				pRef_seq->countRefIndex(alignment,read_len);
				if(alignment.ref_index!=-1){
					output.writeNegSecMappedAlignment(total,alignment);
					aligned_poses++;
				}
			}
		}
		if(mapped){
			mapped_num++;
			continue;
		}
		alignments_size=pSub_read->alignRead(read,pmarks,temp_marks_len,alignments);
		for(i=0;i!=alignments_size;i++){
			alignment.start_pos=alignments[i].ref_pos;
			alignment.cigar=alignments[i].cigar;
			alignment.edit=alignments[i].edit;
			pRef_seq->countRefIndex(alignment,read_len);
			if(alignment.ref_index!=-1){
				output.writePosMappedAlignment(total,read,alignment);
				mapped=true;
				aligned_poses++;
				i++;
				break;
			}
		}
		for(;i!=alignments_size;i++){
			alignment.start_pos=alignments[i].ref_pos;
			alignment.cigar=alignments[i].cigar;
			alignment.edit=alignments[i].edit;
			pRef_seq->countRefIndex(alignment,read_len);
			if(alignment.ref_index!=-1){
				output.writePosSecMappedAlignment(total,alignment);
				aligned_poses++;
			}
		}
		if(mapped){
			mapped_num++;
			continue;
		}
		alignments_size=pSub_read->alignRead(rev_read,nmarks,temp_marks_len,alignments);
		for(i=0;i!=alignments_size;i++){
			alignment.start_pos=alignments[i].ref_pos;
			alignment.cigar=alignments[i].cigar;
			alignment.edit=alignments[i].edit;
			pRef_seq->countRefIndex(alignment,read_len);
			if(alignment.ref_index!=-1){
				output.writeNegMappedAlignment(total,rev_read,alignment);
				mapped=true;
				aligned_poses++;
				i++;
				break;
			}
		}
		for(;i!=alignments_size;i++){
			alignment.start_pos=alignments[i].ref_pos;
			alignment.cigar=alignments[i].cigar;
			alignment.edit=alignments[i].edit;
			pRef_seq->countRefIndex(alignment,read_len);
			if(alignment.ref_index!=-1){
				output.writeNegSecMappedAlignment(total,alignment);
				aligned_poses++;
			}
		}
		if(mapped){
			mapped_num++;
			continue;
		}
		output.writeUnmappedAlignment(total,read);
	}
	input.close();
	output.writeClose();
	cout.precision(2);
	cout<<fixed;
	cout<<blank_line<<endl;
	cout<<"total reads::"<<total<<",mapped reads:"<<mapped_num;
	cout<<",mapped positions:"<<aligned_poses<<endl;
	cout<<"mapped ratio:"<<((float)mapped_num/total)*100<<"%,average positions:";
	if(mapped_num)
		cout<<(float)aligned_poses/mapped_num<<endl;
	else
		cout<<"0"<<endl;
}
double cputime()
{
	struct rusage r;
	getrusage(RUSAGE_SELF, &r);
	return r.ru_utime.tv_sec + r.ru_stime.tv_sec + 1e-6 * (r.ru_utime.tv_usec + r.ru_stime.tv_usec);
}
double realtime()
{
	struct timeval tp;
	struct timezone tzp;
	gettimeofday(&tp, &tzp);
	return tp.tv_sec + tp.tv_usec * 1e-6;
}
void printMarks(bool marks[],int len){
	int i;
	for(i=0;i!=len;i++)
		if(marks[i])
			cout<<'1';
		else
			cout<<'0';
	cout<<endl;
}
