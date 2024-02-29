#include <iostream>
#include <fstream>
#include <vector>
#include "coordinate.h"
#include <algorithm>
#include <string>


using namespace std;

ostream &operator<<(ostream &os, TrainingExample &te)
{
    int n = te.getFeatures().size();
    os << "(";
    for (int i = 0; i < n; i++)
        os << te.getFeature(i) << (i != n - 1 ? ", " : "");
    os << ")";
    return os;
}

string removeSpaces(string word)
{
    string newWord;
    for (int i = 0; i < word.length(); i++)
    {
        if (word[i] != ' ' && word[i] != '\n' && word[i] != '\t')
        {
            newWord += word[i];
        }
    }

    return newWord;
}

int main()
{
    int M = 0, N = 0;
    setprecision(10);
    vector<TrainingExample> ts;
    vector<float> ReplayBuffer;
    ifstream f;
    f.open("doc.txt");
    if (!f.is_open())
    {
        cout << "File not opened" << endl;
        return 1;
    }
    M = 1175;
    N = 8637;
    cout << "M = " << M << ", N = " << N << endl;

    for (std::string s; getline(f, s, ' ');)
    {
      ReplayBuffer.push_back(stof(s));
    }

    for (int i = 0; i < M; i++)
    {
        vector<float> feat(N+1, 1);
        int y = 0;
        for (int j = 0; j < N; j++) 
        {
            feat[j+1] = ReplayBuffer[i*N + j + i];
        }

        y = ReplayBuffer[(i+1) * N];

        TrainingExample te(feat, y);
        ts.push_back(te);
    }
    f.close();

    for (unsigned i = 0; i < ts.size(); i++)
    {
        cout << "Example " << i << ": " << ts[i] << endl;
        cout << "Example " << i << ": " << ts[i].getTarget() << endl;
    }


    SVM hyp(ts);
    vector<float> theta = hyp.coodinateDescent();

    cout << endl;
    for (size_t i = 0; i < theta.size(); i++)
        cout << "th" << i << " = " << theta[i] << " ";
    cout << endl;

    int ex = 0;
    for (size_t i = 0; i < ts.size(); i++)
    {
        if (hyp.predict(ts[i].getFeatures()) - ts[i].getTarget())
        {
            ex += 1;
        }
    }

    cout << ex / ts.size() << endl;

    cout << "Input x1 x2" << endl;
    int x1, x2;
    cin >> x1 >> x2;
    cout << "H = " << (theta[0]+theta[1]*x1+theta[2]*x2) << endl;

    return 0;
}
