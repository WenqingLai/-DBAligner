/*
 * Configure.cpp
 *
 *  Created on: Mar 13, 2014
 *      Author: laiwenqi
 */
#include "Configure.h"
#include<fstream>
#include<iostream>
#include<cstdlib>
#include<iomanip>
#include<sstream>
int FileFormat::fa=0;
int FileFormat::fq=1;
int FileFormat::unknown=2;

string Configure::config_file="./configure";
string Configure::dir="/home/laiwenqi/Documents/Test_Data/";
string Configure::refs_file="Sakai.fa";
string Configure::reads_file="left_reads.fq right_reads.fq";
string Configure::left_reads="";
string Configure::right_reads="";
string Configure::align_file="Aligner1.sam";
int Configure::seed_size=17;
int Configure::max_refs=1000;
int Configure::max_read_len=1000;
float Configure::load_factor=0.75;
int Configure::seed_max_poses=20;
float Configure::global_edit_ratio=0.15;
float Configure::global_gap_ratio=0.08;
float Configure::local_edit_ratio=0.35;
int Configure::local_max_gap=2;
int Configure::local_min_edit=3;
int Configure::max_local_results=200;
int Configure::global_max_edit=0;
int Configure::global_max_gap=0;
float Configure::vote_cutoff=0.2;
bool Configure::report_all=false;
int Configure::reads_fmt=FileFormat::fq;
int Configure::ref_fmt=FileFormat::fa;
bool Configure::is_pe=false;
int Configure::insert_estimate_cutoff=0;
int Configure::forward_chimera_pos=0;
int Configure::backward_chimera_pos=0;
string Configure::read1="";
string Configure::read2="";
string Configure::rev_read1="";
string Configure::rev_read2="";
float Configure::median_factor=0.5;
float Configure::deviation_factor=0;
int Configure::low_bound=-1;
int Configure::high_bound=-1;
bool Configure::bound_designated=true;

bool Configure::initialize(){
	ifstream input;
	string temp;
	const char * temp_chars;
	int i,len,line,begin,end,count;
	bool in_quote,in_digit,in_comment;
	input.open(config_file.c_str());
	if(!input){
		input.clear();
		input.open(("./"+config_file).c_str());
		if(!input){
			cerr<<"Can't find the file 'configure'!"<<endl;
			return false;
		}
	}
	line=0;
	while(getline(input,temp)){
		len=temp.length();
		temp_chars=temp.c_str();
		in_comment=false;
		for(i=0;i!=len;i++){
			if(temp_chars[i]=='#'){
				in_comment=true;
				break;
			}
		}
		if(in_comment)
			continue;
		switch(line){
		case 0:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::dir=temp.substr(begin,end-begin);
			break;
		case 1:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::refs_file=temp.substr(begin,end-begin);
			Configure::ref_fmt=findReferenceFormat(refs_file);
			if(Configure::ref_fmt==FileFormat::unknown){
				cerr<<"can't recognize reference file format.";
				return false;
			}
			break;
		case 2:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::reads_file=temp.substr(begin,end-begin);
			len=reads_file.length();
			for(i=0;i!=len;i++){
				if(isspace(reads_file[i]))
					break;
			}
			if(i==len){
				Configure::is_pe=false;
				Configure::reads_fmt=findReadsFormat(reads_file);
				if(Configure::reads_fmt==FileFormat::unknown){
					cerr<<"can't recognize reads file format.";
					return false;
				}
				break;
			}
			Configure::is_pe=true;
			Configure::left_reads=reads_file.substr(0,i);
			for(;i!=len;i++)
				if(!isspace(reads_file[i]))
					break;
			Configure::right_reads=reads_file.substr(i);
			Configure::reads_fmt=findReadsFormat(left_reads);
			if(Configure::reads_fmt==FileFormat::unknown){
				cerr<<"can't recognize reads file format.";
				return false;
			}
			break;
		case 3:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::align_file=temp.substr(begin,end-begin);
			break;
		case 4:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
				}
			Configure::seed_size=atoi(temp.substr(begin,count).c_str());
			break;
		case 5:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
				}
			Configure::global_edit_ratio=atof(temp.substr(begin,count).c_str());
			break;
		case 6:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::global_gap_ratio=atof(temp.substr(begin,count).c_str());
			break;
		case 7:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::local_edit_ratio=atof(temp.substr(begin,count).c_str());
			break;
		case 8:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
				}
			Configure::local_min_edit=atoi(temp.substr(begin,count).c_str());
			break;
		case 9:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::local_max_gap=atoi(temp.substr(begin,count).c_str());
			break;
		case 10:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::max_local_results=atoi(temp.substr(begin,count).c_str());
			break;
		case 11:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::vote_cutoff=atof(temp.substr(begin,count).c_str());
			break;
		case 12:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			if(atoi(temp.substr(begin,count).c_str()))
				Configure::report_all=true;
			break;
		case 13:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::max_refs=atoi(temp.substr(begin,count).c_str());
			break;
		case 14:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::max_read_len=atoi(temp.substr(begin,count).c_str());
			break;
		case 15:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::insert_estimate_cutoff=atoi(temp.substr(begin,count).c_str());
			break;
		case 16:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::median_factor=atof(temp.substr(begin,count).c_str());
			break;
		case 17:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::deviation_factor=atof(temp.substr(begin,count).c_str());
			break;
		case 18:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			if(atoi(temp.substr(begin,count).c_str())){
				bound_designated=true;
			}
			else{
				bound_designated=false;
				return true;
			}
			break;
		case 19:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::low_bound=atoi(temp.substr(begin,count).c_str());
			break;
		case 20:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::high_bound=atoi(temp.substr(begin,count).c_str());
			break;
		default:
			break;
		}
		line++;
	}
	input.close();
	return true;
}
bool Configure::initialize(const string & file){
	ifstream input;
	string temp;
	const char * temp_chars;
	int i,len,line,begin,end,count;
	bool in_quote,in_digit,in_comment;
	temp=Configure::dir+file;
	input.open(temp.c_str());
	if(!input){
		input.clear();
		input.open(("./"+file).c_str());
		if(!input){
			cerr<<"Can't find the configure file!"<<endl;
			return false;
		}
	}
	line=0;
	while(getline(input,temp)){
		len=temp.length();
		temp_chars=temp.c_str();
		in_comment=false;
		for(i=0;i!=len;i++){
			if(temp_chars[i]=='#'){
				in_comment=true;
				break;
			}
		}
		if(in_comment)
			continue;
		switch(line){
		case 0:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::dir=temp.substr(begin,end-begin);
			break;
		case 1:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::refs_file=temp.substr(begin,end-begin);
			Configure::ref_fmt=findReferenceFormat(refs_file);
			if(Configure::ref_fmt==FileFormat::unknown){
				cerr<<"can't recognize reference file format.";
				return false;
			}
			break;
		case 2:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::reads_file=temp.substr(begin,end-begin);
			len=reads_file.length();
			for(i=0;i!=len;i++){
				if(isspace(reads_file[i]))
					break;
			}
			if(i==len){
				Configure::is_pe=false;
				Configure::reads_fmt=findReadsFormat(reads_file);
				if(Configure::reads_fmt==FileFormat::unknown){
					cerr<<"can't recognize reads file format.";
					return false;
				}
				break;
			}
			Configure::is_pe=true;
			Configure::left_reads=reads_file.substr(0,i);
			for(;i!=len;i++)
				if(!isspace(reads_file[i]))
					break;
			Configure::right_reads=reads_file.substr(i);
			Configure::reads_fmt=findReadsFormat(left_reads);
			if(Configure::reads_fmt==FileFormat::unknown){
				cerr<<"can't recognize reads file format.";
				return false;
			}
			break;
		case 3:
			in_quote=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]=='"'){
					if(!in_quote){
						in_quote=true;
						begin=i+1;
					}
					else{
						end=i;
						break;
					}
				}
			}
			Configure::align_file=temp.substr(begin,end-begin);
			break;
		case 4:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
				}
			Configure::seed_size=atoi(temp.substr(begin,count).c_str());
			break;
		case 5:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
				}
			Configure::global_edit_ratio=atof(temp.substr(begin,count).c_str());
			break;
		case 6:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::global_gap_ratio=atof(temp.substr(begin,count).c_str());
			break;
		case 7:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::local_edit_ratio=atof(temp.substr(begin,count).c_str());
			break;
		case 8:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
				}
			Configure::local_min_edit=atoi(temp.substr(begin,count).c_str());
			break;
		case 9:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::local_max_gap=atoi(temp.substr(begin,count).c_str());
			break;
		case 10:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::max_local_results=atoi(temp.substr(begin,count).c_str());
			break;
		case 11:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::vote_cutoff=atof(temp.substr(begin,count).c_str());
			break;
		case 12:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			if(atoi(temp.substr(begin,count).c_str()))
				Configure::report_all=true;
			break;
		case 13:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::max_refs=atoi(temp.substr(begin,count).c_str());
			break;
		case 14:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::max_read_len=atoi(temp.substr(begin,count).c_str());
			break;
		case 15:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::insert_estimate_cutoff=atoi(temp.substr(begin,count).c_str());
			break;
		case 16:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::median_factor=atof(temp.substr(begin,count).c_str());
			break;
		case 17:
			in_digit=false;
			for(i=0;i!=len;i++){
				if((temp_chars[i]>='0'&&temp_chars[i]<='9')||temp_chars[i]=='.'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
					}
			}
			Configure::deviation_factor=atof(temp.substr(begin,count).c_str());
			break;
		case 18:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			if(atoi(temp.substr(begin,count).c_str())){
				bound_designated=true;
			}
			else{
				bound_designated=false;
				return true;
			}
			break;
		case 19:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::low_bound=atoi(temp.substr(begin,count).c_str());
			break;
		case 20:
			in_digit=false;
			for(i=0;i!=len;i++){
				if(temp_chars[i]>='0'&&temp_chars[i]<='9'){
					if(!in_digit){
						in_digit=true;
						begin=i;
						count=1;
					}
					else
						count++;
				}
			}
			Configure::high_bound=atoi(temp.substr(begin,count).c_str());
			break;
		default:
			break;
		}
		line++;
	}
	input.close();
	return true;
}
int Configure::findReferenceFormat(const string & ref_file){
	int i;
	bool has_dot;
	string temp;
	has_dot=false;
	for(i=ref_file.length()-1;i!=-1;i--){
		if(ref_file[i]=='.'){
			has_dot=true;
			break;
		}
	}
	if(!has_dot)
		return FileFormat::unknown;
	temp=ref_file.substr(i+1);
	if(temp=="fa")
		return FileFormat::fa;
	else if(temp=="fasta")
		return FileFormat::fa;
	else
		return FileFormat::unknown;
}
int Configure::findReadsFormat(const string & reads_file){
	int i;
	bool has_dot;
	string temp;
	has_dot=false;
	for(i=reads_file.length()-1;i!=-1;i--){
		if(reads_file[i]=='.'){
			has_dot=true;
			break;
		}
	}
	if(!has_dot)
		return FileFormat::unknown;
	temp=reads_file.substr(i+1);
	if(temp=="fq")
		return FileFormat::fq;
	else if(temp=="fastq")
		return FileFormat::fq;
	else if(temp=="fa")
		return FileFormat::fa;
	else if(temp=="fasta")
		return FileFormat::fa;
	else
		return FileFormat::unknown;
}
void Configure::printSingleConfiguration(){
	ostringstream temp_out;
	cout<<refs_file;
	if(Configure::is_pe){
		cout<<","<<left_reads<<","<<right_reads;
	}
	else{
		cout<<","<<reads_file;
	}
	cout<<","<<align_file<<endl;
	temp_out<<"seed_size:"<<seed_size;
	cout<<setw(25)<<left<<temp_out.str();
	temp_out.str("");
	temp_out<<"max_read_len:"<<max_read_len;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	temp_out.str("");
	temp_out<<"global_edits_ratio:"<<global_edit_ratio;
	cout<<setw(25)<<left<<temp_out.str();
	temp_out.str("");
	temp_out<<"global_gap_bases_ratio:"<<global_gap_ratio;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	temp_out.str("");
	temp_out<<"local_edits_ratio:"<<local_edit_ratio;
	cout<<setw(25)<<left<<temp_out.str();
	temp_out.str("");
	temp_out<<"max_local_gap_bases:"<<local_max_gap;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	cout<<setw(25)<<left;
	if(Configure::report_all)
		cout<<"report_all:true";
	else
		cout<<"report_all:false";
	temp_out.str("");
	temp_out<<"vote_cutoff:"<<vote_cutoff;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	cout.flush();
}
void Configure::printPEConfiguration(){
	//cout<<dir<<endl;
	ostringstream temp_out;
	cout<<refs_file;
	if(Configure::is_pe){
		cout<<","<<left_reads<<","<<right_reads;
	}
	else{
		cout<<","<<reads_file;
	}
	cout<<","<<align_file<<endl;
	temp_out<<"seed_size:"<<seed_size;
	cout<<setw(25)<<left<<temp_out.str();
	temp_out.str("");
	temp_out<<"max_read_len:"<<max_read_len;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	temp_out.str("");
	temp_out<<"global_edits_ratio:"<<global_edit_ratio;
	cout<<setw(25)<<left<<temp_out.str();
	temp_out.str("");
	temp_out<<"global_gap_bases_ratio:"<<global_gap_ratio;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	temp_out.str("");
	temp_out<<"local_edits_ratio:"<<local_edit_ratio;
	cout<<setw(25)<<left<<temp_out.str();
	temp_out.str("");
	temp_out<<"max_local_gap_bases:"<<local_max_gap;
	cout<<setw(25)<<left<<temp_out.str()<<endl;
	temp_out.str("");
	temp_out<<"vote_cutoff:"<<vote_cutoff;
	cout<<setw(25)<<left<<temp_out.str();
	if(!bound_designated){
		temp_out.str("");
		temp_out<<"median_factor:"<<median_factor;
		cout<<setw(25)<<left<<temp_out.str()<<endl;
		temp_out.str("");
		temp_out<<"lest_stone_pes:"<<insert_estimate_cutoff;
		cout<<setw(25)<<left<<temp_out.str();
		temp_out.str("");
		temp_out<<"deviation_factor:"<<deviation_factor;
		cout<<setw(25)<<left<<temp_out.str()<<endl;
	}
	else{
		cout<<setw(25)<<left<<"designate_insert_scope:true"<<endl;
		temp_out.str("");
		temp_out<<"insert_low_bound:"<<Configure::low_bound;
		cout<<setw(25)<<left<<temp_out.str();
		temp_out.str("");
		temp_out<<"insert_high_bound:"<<Configure::high_bound;
		cout<<setw(25)<<left<<temp_out.str()<<endl;
	}
	cout.flush();
}
