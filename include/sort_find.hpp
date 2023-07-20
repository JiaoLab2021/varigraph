#ifndef sort_find_hpp
#define sort_find_hpp
#include<iostream>
#include<algorithm>
#include <vector>
#include<regex>

using namespace std;

template<typename T>

// ��������
int findPosArray(T ar[], int n, T element)//����Ԫ�ز�����λ���±꣬find(���飬���ȣ�Ԫ��)
{
	int i = 0;
	int index=-1;//ԭʼ�±꣬û�ҵ�Ԫ�ط���-1
	for (i = 0; i <n; i++)
	{
		if (element ==ar[i])
		{
			index=i;//��¼Ԫ���±�
		}
	}
	return index;//�����±�
}

// ��������
int findPosVector(vector <int> input , int number)
{
    vector<int>::iterator iter=std::find(input.begin(),input.end(),number);//���ص���һ��������ָ��
    if(iter == input.end())
    {
        return -1;
    } else
	{
        return std::distance(input.begin(),iter);
    }
}

// ����
bool Reverse(int a,int b)
{
    return a > b; //�������У������Ϊreturn a>b����Ϊ����
}

// ���ֲ���
int search_Binary_left(vector<int>v, int value) // search_Binary_left(����, Ҫ�ҵ���)
{
	int low = 0;
	int high = v.size() - 1;
	int mid = (low + high) / 2;
	while (low <= high)
	{

		if (v[mid] == value)
		{
			return mid;
		}
		else if (value < v[mid])
		{
			high = mid-1;
		}
		else
		{
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}
	if (high < 0)
	{
		high = 0;
	}
	
	return high;
}

int search_Binary_right(vector<int>v, int value) // search_Binary_right(����, Ҫ�ҵ���)
{
	int low = 0;
	int high = v.size() - 1;
	int mid = (low + high) / 2;
	while (low <= high)
	{

		if (v[mid] == value)
		{
			return mid;
		}
		else if (value < v[mid])
		{
			high = mid-1;
		}
		else
		{
			low = mid + 1;
		}
		mid= (low + high) / 2;
	}
	if (low > (v.size() - 1))
	{
		low = (v.size() - 1);
	}
	
	return low;
}

// �ַ����Ķ��ֲ���
int search_Binary_string(vector<string> v, long long int begin, long long int end, string s)
{
	long long int result = -1;
	while(begin <= end)
	{
		int mid = begin+(end-begin)/2;
		if(v[mid] > s)
			end = mid - 1;
		else if(v[mid] < s)
			begin = mid + 1;
		else
		{
			if(result < mid)
				result = mid;
			begin = mid + 1;
		}
	}
	return result;
}

#endif