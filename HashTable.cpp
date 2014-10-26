/*
 * HashTable.cpp
 *
 *  Created on: Apr 21, 2014
 *      Author: laiwenqi
 */

#include "HashTable.h"

HashTable::HashTable(Reference * pRef_seq) {
	ref_str=pRef_seq->ref_str.c_str();
	ref_size=pRef_seq->ref_str.size();
	table_size=ref_size/Configure::load_factor;
	pTable=new int[table_size];
	fill_n(pTable,table_size,-1);
	buckets_capacity=ref_size-Configure::seed_size+1;
	buckets=new Element[buckets_capacity];
	buckets_size=0;
	repeats=new int[buckets_capacity];
	fill_n(repeats,buckets_capacity,-1);
	result_capacity=1<<18;
	presult=new FinalElement[result_capacity];
}

HashTable::~HashTable() {
	// TODO Auto-generated destructor stub
}
int HashTable::hashCode(int begin){
	int i,end;
	register unsigned int hash;
	register unsigned char *p;
	end=begin+Configure::seed_size;
	p=(unsigned char*)(ref_str+begin);
	hash=0;
	for(i=begin;i!=end;i++){
	    hash = 31 * hash + *p;
	    p++;
	}
	return (hash & 0x7FFFFFFF);
}
int HashTable::hashCode(const string & str){
	int i,size;
	register unsigned int hash;
	register unsigned char *p;
	size=str.size();
	hash=0;
	p=(unsigned char*)str.c_str();
	for(i=0;i!=size;i++){
		hash = 31 * hash + *p;
		p++;
	}
	return (hash & 0x7FFFFFFF);
}
bool HashTable::refEqual(int ref_pos1,int ref_pos2){
	int i,j,k;
	for(k=0,i=ref_pos1,j=ref_pos2;k!=Configure::seed_size;k++,i++,j++)
		if(ref_str[i]!=ref_str[j])
			return false;
	return true;
}
bool HashTable::refEqual2(const string & str,int ref_pos){
	int i,k;
	for(i=ref_pos,k=0;k!=Configure::seed_size;k++,i++)
		if(ref_str[i]!=str[k])
			return false;
	return true;
}
bool HashTable::tryToAdd(int ref_pos,int vertex_loc){
	string str;
	int hash_pos,temp_ref_pos,buckets_size_bucket;
	bool found;
	Element element;
	hash_pos=hashCode(ref_pos)%table_size;
	if(pTable[hash_pos]==-1){
		element.ref_pos=ref_pos;
		element.vertex_loc=vertex_loc;
		element.next=buckets_size;
		buckets[buckets_size]=element;
		pTable[hash_pos]=buckets_size;
		buckets_size++;
		return true;
	}
	buckets_size_bucket=buckets[pTable[hash_pos]].next;
	found=false;
	temp_ref_pos=-1;
	while(buckets_size_bucket!=pTable[hash_pos]){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual(temp_ref_pos,ref_pos)){
			found=true;
			break;
		}
		buckets_size_bucket=buckets[buckets_size_bucket].next;
	}
	if(!found){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual(temp_ref_pos,ref_pos))
			found=true;
	}
	if(found){
		repeats[ref_pos]=temp_ref_pos;
		buckets[buckets_size_bucket].ref_pos=ref_pos;
	}
	else{
		element.ref_pos=ref_pos;
		element.vertex_loc=vertex_loc;
		element.next=buckets[pTable[hash_pos]].next;
		buckets[buckets_size]=element;
		buckets[pTable[hash_pos]].next=buckets_size;
		pTable[hash_pos]=buckets_size;
		buckets_size++;
	}
	return !found;
}
int HashTable::get(int ref_pos){
	int hash_pos,temp_ref_pos,buckets_size_bucket;
	hash_pos=hashCode(ref_pos)%table_size;
	buckets_size_bucket=buckets[pTable[hash_pos]].next;
	while(buckets_size_bucket!=pTable[hash_pos]){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual(temp_ref_pos,ref_pos))
			break;
		buckets_size_bucket=buckets[buckets_size_bucket].next;
		}
	return buckets[buckets_size_bucket].vertex_loc;
}
int HashTable::get(const string & str){
	int hash_pos,temp_ref_pos,buckets_size_bucket;
	hash_pos=hashCode(str)%table_size;
	buckets_size_bucket=buckets[pTable[hash_pos]].next;
	while(buckets_size_bucket!=pTable[hash_pos]){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual2(str,temp_ref_pos))
				break;
		buckets_size_bucket=buckets[buckets_size_bucket].next;
	}
	return buckets[buckets_size_bucket].vertex_loc;
}
int HashTable::getRepeatHead(const string & str){
	int hash_pos,temp_ref_pos,buckets_size_bucket;
	bool found;
	hash_pos=hashCode(str)%table_size;
	if(pTable[hash_pos]==-1)
		return -1;
	buckets_size_bucket=buckets[pTable[hash_pos]].next;
	found=false;
	temp_ref_pos=-1;
	while(buckets_size_bucket!=pTable[hash_pos]){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual2(str,temp_ref_pos)){
			found=true;
			break;
		}
		buckets_size_bucket=buckets[buckets_size_bucket].next;
	}
	if(!found){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual2(str,temp_ref_pos))
			found=true;
	}
	if(found)
		return temp_ref_pos;
	return -1;
}
int HashTable::verifyResult(const string & str,int poses[]){
	int head,repeat_head,k,i,times,start,pos,ref_pos,remain;
	bool is_novel;
	int buckets_size_pos,ref_read_size;
	ref_read_size=str.size();
	repeat_head=getRepeatHead(str.substr(0,Configure::seed_size));
	if(repeat_head==-1)
		return 0;
	k=repeat_head;
	i=0;
	while(k!=-1){
		presult[i].ref_pos=k+Configure::seed_size;
		presult[i].next=i+1;
		i++;
		if(i==result_capacity)
			break;
		k=repeats[k];
	}
	presult[i-1].next=-1;
	head=0;
	times=str.length()/Configure::seed_size;
	start=0;
	for(i=1;i!=times;i++){
		start+=Configure::seed_size;
		repeat_head=getRepeatHead(str.substr(start,Configure::seed_size));
		k=repeat_head;
		pos=head;
		is_novel=true;
		while(pos!=-1&&k!=-1){
			ref_pos=presult[pos].ref_pos;
			if(ref_pos>k)
				k=repeats[k];
			else if(ref_pos==k){
				k=repeats[k];
				if(is_novel){
					is_novel=false;
					head=pos;
					repeat_head=pos;
				}
				else{
				    presult[repeat_head].next=pos;
				    repeat_head=pos;
				}
				presult[pos].ref_pos+=Configure::seed_size;
				pos=presult[pos].next;
			}
			else
				pos=presult[pos].next;
		}
		if(is_novel)
			return 0;
		presult[repeat_head].next=-1;
	}
	remain=Configure::seed_size-str.length()%Configure::seed_size;
	if(remain==Configure::seed_size){
		pos=head;
		buckets_size_pos=0;
		while(pos!=-1){
			poses[buckets_size_pos++]=presult[pos].ref_pos-ref_read_size;
			if(buckets_size_pos==Configure::max_local_results)
				break;
			pos=presult[pos].next;
		}
		return buckets_size_pos;
	}
	start=str.length()-Configure::seed_size;
	repeat_head=getRepeatHead(str.substr(start,Configure::seed_size));
	k=repeat_head;
	pos=head;
	is_novel=true;
	while(pos!=-1&&k!=-1){
		ref_pos=presult[pos].ref_pos;
		ref_pos-=remain;
		if(ref_pos>k)
			k=repeats[k];
		else if(ref_pos==k){
			k=repeats[k];
			if(is_novel){
				is_novel=false;
				head=pos;
				repeat_head=pos;
			}
			else{
				presult[repeat_head].next=pos;
				repeat_head=pos;
			}
			presult[pos].ref_pos=ref_pos+Configure::seed_size;
			pos=presult[pos].next;
		}
		else
			pos=presult[pos].next;
	}
	if(is_novel)
		return 0;
	presult[repeat_head].next=-1;
	pos=head;
	buckets_size_pos=0;
	while(pos!=-1){
		poses[buckets_size_pos++]=presult[pos].ref_pos-ref_read_size;
		if(buckets_size_pos==Configure::max_local_results)
			break;
		pos=presult[pos].next;
	}
	return buckets_size_pos;
}
int HashTable::lookupString(const string & str,int poses[],int max_size){
	int head,repeat_head,k,i,times,start,pos,ref_pos,remain;
	bool is_novel;
	int buckets_size_pos,ref_read_size;
	buckets_size_pos=0;
	ref_read_size=str.size();
	repeat_head=getRepeatHead(str.substr(0,Configure::seed_size));
	if(repeat_head==-1)
		return 0;
	k=repeat_head;
	i=0;
	while(k!=-1){
		presult[i].ref_pos=k+Configure::seed_size;
		presult[i].next=i+1;
		i++;
		if(i==result_capacity)
			break;
		k=repeats[k];
	}
	presult[i-1].next=-1;
	head=0;
	times=str.length()/Configure::seed_size;
	start=0;
	for(i=1;i!=times;i++){
		start+=Configure::seed_size;
		repeat_head=getRepeatHead(str.substr(start,Configure::seed_size));
		k=repeat_head;
		pos=head;
		is_novel=true;
		while(pos!=-1&&k!=-1){
			ref_pos=presult[pos].ref_pos;
			if(ref_pos>k)
				k=repeats[k];
			else if(ref_pos==k){
				k=repeats[k];
				if(is_novel){
					is_novel=false;
					head=pos;
					repeat_head=pos;
				}
				else{
				    presult[repeat_head].next=pos;
				    repeat_head=pos;
				}
				presult[pos].ref_pos+=Configure::seed_size;
				pos=presult[pos].next;
			}
			else
				pos=presult[pos].next;
		}
		if(is_novel)
			return 0;
		presult[repeat_head].next=-1;
	}
	remain=Configure::seed_size-str.length()%Configure::seed_size;
	if(remain==Configure::seed_size){
		pos=head;
		while(pos!=-1){
			poses[buckets_size_pos++]=presult[pos].ref_pos-ref_read_size;
			if(buckets_size_pos==max_size)
				break;
			pos=presult[pos].next;
		}
		return buckets_size_pos;
	}
	start=str.length()-Configure::seed_size;
	repeat_head=getRepeatHead(str.substr(start,Configure::seed_size));
	k=repeat_head;
	pos=head;
	is_novel=true;
	while(pos!=-1&&k!=-1){
		ref_pos=presult[pos].ref_pos;
		ref_pos-=remain;
		if(ref_pos>k)
			k=repeats[k];
		else if(ref_pos==k){
			k=repeats[k];
			if(is_novel){
				is_novel=false;
				head=pos;
				repeat_head=pos;
			}
			else{
				presult[repeat_head].next=pos;
				repeat_head=pos;
			}
			presult[pos].ref_pos=ref_pos+Configure::seed_size;
			pos=presult[pos].next;
		}
		else
			pos=presult[pos].next;
	}
	if(is_novel)
		return 0;
	presult[repeat_head].next=-1;
	pos=head;
	while(pos!=-1){
		poses[buckets_size_pos++]=presult[pos].ref_pos-ref_read_size;
		if(buckets_size_pos==max_size)
			break;
		pos=presult[pos].next;
	}
	return buckets_size_pos;
}
bool HashTable::containsKey(const string & str){
	int hash_pos,buckets_size_bucket,temp_ref_pos;
	bool found;
	hash_pos=hashCode(str)%table_size;
	if(pTable[hash_pos]==-1)
		return false;
	buckets_size_bucket=buckets[pTable[hash_pos]].next;
	found=false;
	while(buckets_size_bucket!=pTable[hash_pos]){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual2(str,temp_ref_pos)){
			found=true;
			break;
		}
		buckets_size_bucket=buckets[buckets_size_bucket].next;
	}
	if(!found){
		temp_ref_pos=buckets[buckets_size_bucket].ref_pos;
		if(refEqual2(str,temp_ref_pos))
			found=true;
	}
	return found;
}
