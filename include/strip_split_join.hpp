#ifndef STRIP_SPLIT_JOIN_HPP
#define STRIP_SPLIT_JOIN_HPP
#include <string>
#include <vector>

using namespace std;


/**
 * Remove special characters at both ends of the string
 * str -> String
 * ch -> Special characters
**/
string strip(const string& str, char ch=' ');


/**
 * Split the string
 * str -> String
 * delim -> Split characters
**/
vector<string> split(const string& str, const string & delim);

/**
 * String join merging
 * val -> The container to be joined
 * delim -> Characters added between container elements
**/
template<typename T>
string join(vector<T>& val, string delim);

/**
 * String join merging
 * val -> The container to be joined
 * delim -> Characters added between container elements
**/
template<typename T>
string join(const vector<T>& val, string delim);

#endif