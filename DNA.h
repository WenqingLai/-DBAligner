/*
 * DNA.h
 *
 *  Created on: Mar 9, 2014
 *      Author: laiwenqi
 */

#ifndef DNA_H_
#define DNA_H_
#include<string>
using namespace std;
class DNA {
public:
	DNA();
	virtual ~DNA();
	static void reverseComplete(char chars[],int len);
	static string reverse(const string & str);
	static string reverseComplete(const string & str);
};

#endif /* DNA_H_ */
