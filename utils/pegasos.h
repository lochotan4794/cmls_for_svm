#ifndef GRADIENT_H
#define GRADIENT_H

#include <vector>
#include <iostream>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <stdlib.h> /* abs */
#include <numeric>
#include <vector>
#include <cstdlib> // std::rand, std::srand
#include <ctime>   // std::time
#include <cstdlib> // std::rand, std::srand
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <numeric>
#include <vector>
#include <cmath>

using std::cout;
using std::endl;
using std::max;
using std::swap;
using std::vector;
using std::min;



class TrainingExample
{
private:
    vector<float> features;
    int target;

public:
    TrainingExample(vector<float> &feat, int tar)
    {
        features = feat;
        target = tar;
    }
    float getFeature(int i) { return features.at(i); }
    vector<float> &getFeatures() { return features; }
    int getTarget() { return target; }
};

class SVM
{
private:
    vector<float> theta;
    vector<float> z;
    vector<TrainingExample> ts;
    unsigned mExamples, nFeatures;
    float C;

    float J()
    {
        float sum = 0.0;
        for (unsigned i = 0; i < mExamples; i++)
        {
            float diff = dot(theta, ts[i].getFeatures());
            sum += fmaxf(0.0, 1.0 - diff * ts[i].getTarget());
        }
        return sum / mExamples + C * pow_theta(theta) * 0.5;
    }

    float dot(vector<float> &ntheta, vector<float> &features, float C = 1)
    {
        // cout << "H(x) = ";
        float sum = 0.0;
        for (unsigned i = 0; i < nFeatures; i++)
        {
            // cout << ntheta[i] << "*" << features[i] << " ";
            sum += ntheta[i] * features[i];
        }
        // cout << " = " << sum << endl;
        return sum;
    }

    float pow_theta(vector<float> &ntheta)
    {
        float sum = 0.0;
        for (unsigned i = 0; i < nFeatures; i++)
        {
            // cout << ntheta[i] << "*" << features[i] << " ";
            sum += pow(ntheta[i], 2);
        }
        return sum;
    }

    float I(vector<float> &ntheta, vector<float> &features, int target)
    {
        float sum = 0.0;
        for (unsigned i = 0; i < nFeatures; i++)
        {
            // cout << ntheta[i] << "*" << features[i] << " ";
            sum += ntheta[i] * features[i];
        }

        if (sum * target < 1.0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    vector<float> MUL(float a, vector<float> &b)
    {
        vector<float> result;
        for (unsigned i = 0; i < b.size(); i++)
        {
            result.push_back(b[i] * a);
        }

        return result;
    }

    vector<float> SUM(vector<float> &a, vector<float> &b)
    {
        vector<float> result;
        for (unsigned i = 0; i < b.size(); i++)
        {
            result.push_back(b[i] + a[i]);
        }

        return result;
    }

public:
    SVM(vector<TrainingExample> &examples)
    {
        nFeatures = 0;
        ts = examples;
        mExamples = examples.size();
        if (mExamples > 0)
            nFeatures = examples[0].getFeatures().size();
        theta = vector<float>(nFeatures, 0.0);
    }

    vector<float> fit(int T, int S, float lamda, int k = 100)
    {
        const float alpha = 0.0000001;
        const float eps = 0.00001;
        bool converge = false;
        int debug = 0;
        float diff = 0;
        diff = J();
        vector<float> A;
        C = lamda;

        cout << "J(theta) = " << diff << endl
             << endl;

        // Providing a seed value
        srand((unsigned)time(NULL));

        for (unsigned t = 1; t < T; t++)
        {
            // cout << "Using example" << i << endl;
            for (unsigned it = 1; it < k; it++)
            {
                int i = (rand() % mExamples);

                float mu = 1.0 / (lamda * float(t));
                TrainingExample s = ts[i];

                if (I(theta, s.getFeatures(), s.getTarget()))
                {
                    vector<float> a = MUL((1.0 - mu * lamda), theta);
                    vector<float> b = MUL((mu / k) * s.getTarget(), s.getFeatures());
                    theta = SUM(a, b);
                }
                // else
                // {
                //     theta = MUL((1 - mu * lamda), theta);
                // }
            }

            theta = MUL(fmin(1.0, (1/sqrt(lamda))/sqrt(dot(theta, theta))), theta);
            cout << "t = " << t << endl;
            cout << "J(theta) = " << J() << endl
                 << endl;
        }
        return theta;
    }

    float predict(vector<float> &features)
    {
        // cout << "H(x) = ";
        float sum = 0.0;
        for (unsigned i = 0; i < nFeatures; i++)
        {
            // cout << ntheta[i] << "*" << features[i] << " ";
            sum += theta[i] * features[i];
        }
        // cout << " = " << sum << endl;
        if (sum > 0.0)
        {
            return 1.0;
        }

        return -1.0;
    }
};

#endif
