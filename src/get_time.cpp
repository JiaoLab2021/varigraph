// g++ -c src/get_time.cpp -std=c++17 -O3 -march=native
#include "../include/get_time.hpp"

using namespace std;

string getTime()
{
    time_t timep;
    time (&timep);
    char tmp[64];
    strftime(tmp, sizeof(tmp), "%Y-%m-%d %H:%M:%S",localtime(&timep));
    return tmp;
}