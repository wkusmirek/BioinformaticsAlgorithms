/*
 Petar 'PetarV' Velickovic
 Algorithm: Needleman-Wunsch

 g++ Needleman-Wunsch.cpp -o nw


 valgrind --tool=callgrind ./nw
 kcachegrind

 valgrind --tool=memcheck ./nw

 valgrind --tool=massif --heap=yes --stacks=no ./nw
 massif-visualizer
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <algorithm>
#include <queue>
#include <stack>
#include <set>
#include <map>
#include <complex>

#define DPRINTC(C) printf(#C " = %c\n", (C))
#define DPRINTS(S) printf(#S " = %s\n", (S))
#define DPRINTD(D) printf(#D " = %d\n", (D))
#define DPRINTLLD(LLD) printf(#LLD " = %lld\n", (LLD))
#define DPRINTLF(LF) printf(#LF " = %.5lf\n", (LF))

using namespace std;
typedef long long lld;
typedef unsigned long long llu;

int match_score = 2, mismatch_score = 1, gap_score = 1;
string A_1 = "A";
string B_1 = "B";
string A_100 = "CATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCT";
string B_100 = "ACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCT";
string A_200 = "CATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCT";
string B_200 = "ACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCT";
string A_400 = "ACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCT";
string B_400 = "ACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTACGCTCTCACCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCTCATGTATCCT";
string A = A_400, B = B_400;
const int n = 400, m = 400;
int dp[n][m];

/*
 Needleman-Wunsch algorithm for determining the optimal alignment between two strings
 assuming a given score for hits, gaps and mismatches.
 Complexity: O(n * m) time, O(n * m) memory
*/

inline int needleman_wunsch()
{
    for (int i=0;i<=n;i++) dp[i][0] = dp[0][i] = -i * gap_score;
    for (int i=1;i<=n;i++)
    {
        for (int j=1;j<=m;j++)
        {
            int S = (A[i-1] == B[j-1]) ? match_score : -mismatch_score;
            dp[i][j] = max(dp[i-1][j-1] + S, max(dp[i-1][j] - gap_score, dp[i][j-1] - gap_score));
        }
    }
    return dp[n][m];
}

inline pair<string, string> get_optimal_alignment()
{
    string retA, retB;
    stack<char> SA, SB;
    int ii = n, jj = m;
    while (ii != 0 || jj != 0)
    {
        if (ii == 0)
        {
            SA.push('-');
            SB.push(B[jj-1]);
            jj--;
        }
        else if (jj == 0)
        {
            SA.push(A[ii-1]);
            SB.push('-');
            ii--;
        }
        else
        {
            int S = (A[ii-1] == B[jj-1]) ? match_score : -mismatch_score;
            if (dp[ii][jj] == dp[ii-1][jj-1] + S)
            {
                SA.push(A[ii-1]);
                SB.push(B[jj-1]);
                ii--; jj--;
            }
            else if (dp[ii-1][jj] > dp[ii][jj-1])
            {
                SA.push(A[ii-1]);
                SB.push('-');
                ii--;
            }
            else
            {
                SA.push('-');
                SB.push(B[jj-1]);
                jj--;
            }
        }
    }
    while (!SA.empty())
    {
        retA += SA.top();
        retB += SB.top();
        SA.pop();
        SB.pop();
    }
    return make_pair(retA, retB);
}

int main()
{
    int** arr = (int**)malloc(n * sizeof(int*));
    for (int i = 0; i < n; i++)
        arr[i] = (int*)malloc(m * sizeof(int));
    printf("%d\n",needleman_wunsch());
    pair<string, string> alignment = get_optimal_alignment();
    printf("%s\n%s\n", alignment.first.c_str(), alignment.second.c_str());
    return 0;
}
