#ifndef NeedlemanWunsch_hpp
#define NeedlemanWunsch_hpp
#include <iostream>
#include <iomanip>
#include <string>

using namespace std;

namespace alignment
{
    long long int max_three(long long int a, long long int b, long long int c);
    tuple<string,string,string> alignment(string seq1, string seq2, int match, int mismatch,
                                        int indel);

    // ����ṹ��
    struct Score
    {
        int match;
        int mismatch;
        int indel;
    };

    // // ���ֵ
    // int score_num = 2;
    // Score score[score_num] = 
    // {
    //     {9, -6, -2},
    //     {9, -3, -2}
    // };

    // tie(seq1_ali, seq2_ali, seq_mid) = alignment(seq1, seq2, match, mismatch, indel);
    // cout << endl << seq1_ali << endl << seq_mid << endl << seq2_ali << endl <<  endl << endl;

    long long int max_three(long long int a, long long int b, long long int c)
    {
        long long int max_num;
        if (a<b)
        {
            max_num = b;
        }
        else
        {
            max_num = a;
        }

        if (max_num < c)
        {
            max_num = c;
        }
        
        return max_num;
    }

    tuple<string,string,string> alignment(string seq1, string seq2, int match, int mismatch, int indel)
    {
        // ������ά����
        const long long int matrix_x = seq1.length() + 1;
        const long long int matrix_y = seq2.length() + 1;

        // cout << seq1.length() << " " << seq2.length() << endl;

        // ��ջ�Ͽ��ٶ�ά���� 1M-2M
        // long long int matrix[matrix_x][matrix_y] = {};

        // ���ڴ��Ͽ��ٶ�̬��ά����
        long long int **matrix;
        matrix = new long long int *[matrix_x];
        for (long long int i = 0; i < matrix_x; i++)
        {
            matrix[i] = new long long int[matrix_y];
        }
        
        // long long int (*matrix)[matrix_y] = new long long int[matrix_x][matrix_y];

        for (long long int i = 0; i < matrix_x; i++)
        {
            for (long long int j = 0; j < matrix_y; j++)
            {
                if (i == 0)
                {
                    if (j == 0)
                    {
                        matrix[i][j] = 0;
                    }
                    else
                    {
                        matrix[i][j] = matrix[i][j-1] + indel;
                    }
                    
                }
                else if (j == 0)
                {
                    matrix[i][j] = matrix[i-1][j] + indel;
                }
                else // Ѱ�����·��
                {   
                    long long int ins = matrix[i][j-1] + indel;
                    long long int del = matrix[i-1][j] + indel;
                    long long int match_score;
                    long long int max_num;
                    if (seq2[j-1] == seq1[i-1])
                    {
                        match_score = matrix[i-1][j-1] + match;
                    }
                    else
                    {
                        match_score = matrix[i-1][j-1] + mismatch;
                    }
                    max_num = max_three(ins, del, match_score);
                    matrix[i][j] = max_num; // ������ֵΪ���ֵ
                                
                }
                // // ��ӡ�ؾ���
                // cout << setw(5) << matrix[i][j];
                
            }
            // cout << endl;
        }
    
        // ����
        long long int seq1_len = matrix_x - 1;
        long long int seq2_len = matrix_y - 1;
        string seq1_ali{};
        string seq2_ali{};
        
        for (long long int i = 0; i < max(matrix_x, matrix_y); i++)
        {
            while (seq1_len > 0 || seq2_len > 0)
            {
                // ���ж���û�������Ѿ���������
                if (seq2_len == 0)
                {
                    seq1_ali = seq1[seq1_len-1] + seq1_ali;
                    seq2_ali = "-" + seq2_ali;
                    seq1_len -= 1;
                    continue;
                }
                if (seq1_len == 0)
                {
                    seq1_ali = "-" + seq2_ali;
                    seq2_ali = seq2[seq2_len-1] + seq2_ali;
                    seq2_len -= 1;
                    continue;
                }

                // ������ƥ��
                if (matrix[seq1_len][seq2_len] - 
                    matrix[seq1_len-1][seq2_len-1] == match && 
                    seq1[seq1_len-1] == seq2[seq2_len-1])
                {
                seq1_ali = seq1[seq1_len-1] + seq1_ali;
                seq2_ali = seq2[seq2_len-1] + seq2_ali;
                seq1_len -= 1;
                seq2_len -= 1;
                }
                
                // �����в�ƥ��
                else if (matrix[seq1_len][seq2_len] - 
                        matrix[seq1_len-1][seq2_len-1] == mismatch && 
                        seq1[seq1_len-1] != seq2[seq2_len-1])
                {
                seq1_ali = seq1[seq1_len-1] + seq1_ali;
                seq2_ali = seq2[seq2_len-1] + seq2_ali;
                seq1_len -= 1;
                seq2_len -= 1;
                }
                
                // y���seq��ins
                else if (matrix[seq1_len][seq2_len] - 
                        matrix[seq1_len][seq2_len-1] == indel)
                {
                seq2_ali = seq2[seq2_len-1] + seq2_ali;
                seq1_ali = "-" + seq1_ali;
                seq2_len -= 1;
                }

                // x���seq��ins
                else if (matrix[seq1_len][seq2_len] - 
                        matrix[seq1_len-1][seq2_len] == indel)
                {
                seq1_ali = seq1[seq1_len-1] + seq1_ali;
                seq2_ali = "-" + seq2_ali;
                seq1_len -= 1;
                }
            }
        }
                
        // alignment�м��ֶ�
        string seq_mid{};
        for (long long int i = 0; i < seq1_ali.length(); i++)
        {
            if (seq1_ali[i] == seq2_ali[i])
            {
                seq_mid += "|";
            }
            else
            {
                seq_mid += " ";
            }
        }

        // �ͷŶ�ά�����ڴ�
        for(long long int i = 0; i < matrix_x; i++)
        {
            delete[] matrix[i];
        }
        delete[] matrix;
        
        return make_tuple(seq1_ali, seq2_ali, seq_mid);
    }
}

#endif