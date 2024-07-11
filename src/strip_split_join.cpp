// g++ -c strip_split_join.cpp -std=c++17 -O3 -march=native

#include "../include/strip_split_join.hpp"

using namespace std;


/**
 * Remove special characters at both ends of the string
 * str -> String
 * ch -> Special characters
**/
string strip(const string & str, char ch)
{
	int i = 0;
	while (str[i] == ch)  // The number of ch characters in the header is i
		i++;
	int j = str.size() - 1;
	while (str[j] == ch)  // The number of ch characters at the end is str.size() -1-j
		j--;		
	return str.substr(i, j+1-i);
}


/**
 * Split the string
 * str -> String
 * delim -> Split characters
**/
vector<string> split(const string & str, const string & delim)
{
	vector<string> res;  // Stores the split substring in a vector
	if("" == str) return  res;  // Return an empty vector if it is empty
	
	string strs = str + delim;  // Expands the string to retrieve the last delimited string
	size_t pos;
	size_t size = strs.size();
 
	for (int i = 0; i < size; ++i) {
		pos = strs.find(delim, i);  // pos indicates the position where the delimiter first appears. The string from i to pos is the delimited string
		if( pos < size)  // If yes, if no separator is found, pos is string::npos
		{
			string s = strs.substr(i, pos - i);  // A substring with length pos-i starting from i
			res.push_back(s);  // The string cut between two consecutive Spaces is an empty string, there is no judgment here whether s is empty, so the output of empty characters in the final result,
			i = pos + delim.size() - 1;
		}
	}
	return res;	
}


/**
 * String join merging
 * val -> The container to be joined
 * delim -> Characters added between container elements
**/
template<typename T>
string join(vector<T>& val, string delim) {
    string str;
    int vecSize = val.size();
    int index = 0;
    for (const auto&iter : val) {
        if constexpr (std::is_same_v<T, std::string>) {
            str += iter;
        } else {
            str += to_string(iter);
        }
        if (index != vecSize - 1) {
            str += delim;
        }
        index++;
    }
    return str;
}


/**
 * String join merging
 * val -> The container to be joined
 * delim -> Characters added between container elements
**/
template<typename T>
string join(const vector<T>& val, string delim) {
    string str;
    int vecSize = val.size();
    int index = 0;
    for (const auto&iter : val) {
        if constexpr (std::is_same_v<T, std::string>) {
            str += iter;
        } else {
            str += to_string(iter);
        }
        if (index != vecSize - 1) {
            str += delim;
        }
        index++;
    }
    return str;
}



template string join(vector<string>& val, string delim);
template string join(vector<int>& val, string delim);
template string join(vector<uint16_t>& val, string delim);
template string join(vector<uint32_t>& val, string delim);
template string join(vector<uint64_t>& val, string delim);
template string join(vector<long int>& val, string delim);
template string join(vector<float>& val, string delim);

template string join(const vector<string>& val, string delim);
template string join(const vector<int>& val, string delim);
template string join(const vector<uint16_t>& val, string delim);
template string join(const vector<uint32_t>& val, string delim);
template string join(const vector<uint64_t>& val, string delim);
template string join(const vector<long int>& val, string delim);
template string join(const vector<float>& val, string delim);