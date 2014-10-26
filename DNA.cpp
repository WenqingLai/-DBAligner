/*
 * DNA.cpp
 *
 *  Created on: Mar 9, 2014
 *      Author: laiwenqi
 */

#include "DNA.h"
#include<algorithm>
#include<iostream>
DNA::DNA() {
	// TODO Auto-generated constructor stub

}

DNA::~DNA() {
	// TODO Auto-generated destructor stub
}
void DNA::reverseComplete(char chars[],int len){
	int i,j;
	char c,temp;
	i=0;j=len-1;
	while(i<=j){
		switch(chars[i]){
			case 'A':
				c='T';
				break;
			case 'T':
				c='A';
				break;
			case 'G':
				c='C';
				break;
			case 'C':
				c='G';
				break;
			default:
				c='N';
				break;
			}
		temp=chars[j];
		chars[j]=c;
		switch(temp){
			case 'A':
				c='T';
				break;
			case 'T':
				c='A';
				break;
			case 'G':
				c='C';
				break;
			case 'C':
				c='G';
				break;
			default:
				c='N';
				break;
			}
		chars[i]=c;
		i++;
		j--;
		}
}
string DNA::reverseComplete(const string & str){
	int i,j;
	char c;
	string buffer;
	buffer.resize(str.length());
	for(j=0,i=str.length()-1;i!=-1;i--,j++){
		switch(str[i]){
			case 'A':
				c='T';
				break;
			case 'T':
				c='A';
				break;
			case 'G':
				c='C';
				break;
			case 'C':
				c='G';
				break;
			default:
				c='N';
				break;
		}
		buffer[j]=c;
	}
	return buffer;
}
string DNA::reverse(const string & str){
	int i,j,str_len;
	string temp_str;
	str_len=str.length();
	temp_str.resize(str_len);
	i=0;;
	for(j=0,i=str_len-1;i!=-1;i--,j++)
		temp_str[i]=str[j];
	return temp_str;
}

