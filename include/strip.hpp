#ifndef strip_hpp
#define strip_hpp
#include<string>

using namespace std;

// ���ַ������˵������ַ�ȥ��
string strip(const string &str,char ch=' ')
{
	//��ȥstr���˵�ch�ַ�
	int i = 0;
	while (str[i] == ch)// ͷ��ch�ַ������� i
		i++;
	int j = str.size() - 1;
	while (str[j] == ch ) //
		j--;		
	return str.substr(i, j+1 -i );
}

#endif