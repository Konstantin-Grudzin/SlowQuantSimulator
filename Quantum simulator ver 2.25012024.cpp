#include "qMatrix.h"
#include<ctime>

//class qProgram
//{
//    vector<string> commands;
//    vector<vector<string>> param;
//    unordered_map<string, int> m;
//public:
//    friend istream& operator<<(istream& in, qProgram& P)
//    {
//        string s;
//        while (true)
//        {
//            //---------------------
//            getline(in, s);
//            if (s == "end") break;
//            //---------------------
//
//            int i = 0;
//            string com;
//            vector<string> par;
//            for (; s[i] != ' '; ++i)
//                com += s[i];
//            P.commands.emplace_back(com);
//            i++;
//            string p;
//            for (; i < s.size(); ++i)
//            {
//                if (s[i] == ',')
//                {
//                    par.emplace_back(p);
//                    P.m[p] = -1;
//                    p = "";
//                }
//                else
//                    p += s[i];
//            }
//            par.emplace_back(p);
//            P.param.emplace_back(par);
//        }
//
//        cout << "Paramaters fom file (Y/N)?";
//        char chs;
//        cin >> chs;
//        if (chs == 'Y')
//            P.Init(in);
//        else
//            P.Init(cin);
//    }
//
//    void Init(istream& in)
//    {
//        for (auto& x : m)
//        {
//            cout << x.first << "? ";
//            in >> x.second;
//        }
//    }
//
//};


//Реализация с 3n кубитами
void Three_n_Qubit_Addition()
{
    qMatrix q(10);
    cout << "Input two nubers with space:";

    int a, b; cin >> a >> b;
    for (int i = 0; i < 3; ++i)
    {
        if ((a >> i) & 1)
        {
            q.X(i);
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        if ((b >> i) & 1)
        {
            q.X(6+i);
        }
    }
    cout << "Thats how it's looks like in quants: "; q.OutReal();
    CARRY(q, 3, 0, 6, 4);
    CARRY(q, 4, 1, 7, 5);
    CARRY(q, 5, 2, 8, 9);
    q.CNOT(2, 8);
    SUM(q, 5, 2, 8);
    RCARRY(q, 4, 1, 7, 5);
    SUM(q, 4, 1, 7);
    RCARRY(q, 3, 0, 6, 4);
    SUM(q, 3, 0, 6);
    
    cout << "Summary of our addition: ";
    q.OutReal();
}

//Реализация с 2n кубитами
void Two_n_Qubit_Addition()
{
    qMatrix q(7);
    cout << "Input two nubers with space:";

    int a, b; cin >> a >> b;
    for (int i = 0; i < 3; ++i)
    {
        if ((b >> i) & 1)
        {
            q.X(3+i);
        }
    }

    cout << "Thats how it's looks like in quants: "; q.OutReal();
    CARRY2(q, 3, 0, 6, 4,a);
    CARRY2(q, 4, 1, 7, 5,a);
    CARRY2(q, 5, 2, 8, 9,a);
    if((a>>2)&1)
        q.X(5);
    SUM2(q, 5, 2, 8,a);
    RCARRY2(q, 4, 1, 7, 5,a);
    SUM2(q, 4, 1, 7,a);
    RCARRY2(q, 3, 0, 6, 4,a);
    SUM2(q, 3, 0, 6,a);
    cout << "Summary of our addition: ";
    q.OutReal();
}
    
int main()
{
    
    ios_base::sync_with_stdio(0);
    cout.tie(NULL);
    time_t start, end;
    time(&start);
    qMatrix q(17);
    q[0] = 0;
    q[2] = 1;
    cout<<q.fast_QSHOR(8,7*13,7,1,7,8,15,16)<<"\n";
    time(&end);
    double seconds = difftime(end, start);
    cout << seconds;
}