#ifndef strip_hpp
#define strip_hpp
#include<string>

using namespace std;

// 将字符串两端的特殊字符去掉
string strip(const string &str,char ch=' ')
{
	//除去str两端的ch字符
	int i = 0;
	while (str[i] == ch)// 头部ch字符个数是 i
		i++;
	int j = str.size() - 1;
	while (str[j] == ch ) //
		j--;		
	return str.substr(i, j+1 -i );
}

#endif