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

using std::cout;
using std::endl;
using std::max;
using std::swap;
using std::vector;

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

    float H(vector<float> &ntheta, vector<float> &features, float C = 1)
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

    float Predict(vector<float> &ntheta, vector<float> &features)
    {
        // cout << "H(x) = ";
        float sum = 0.0;
        for (unsigned i = 0; i < nFeatures; i++)
        {
            // cout << ntheta[i] << "*" << features[i] << " ";
            sum += ntheta[i] * features[i];
        }
        // cout << " = " << sum << endl;
        if (sum > 0.0)
        {
            return 1.0;
        }

        return -1.0;
    }

    float J()
    {
        float sum = 0.0;
        for (unsigned i = 0; i < mExamples; i++)
        {
            float diff = H(theta, ts[i].getFeatures());
            sum += pow(fmaxf(0.0, 1 - diff * ts[i].getTarget()), 2);
        }
        return sum / mExamples + C * pow_theta(theta);
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

    float _b_func(vector<float> &ntheta, vector<float> &features, int target)
    {
        // float sum = 0.0;
        // for (unsigned i = 0; i < nFeatures; i++)
        // {
        //     // cout << ntheta[i] << "*" << features[i] << " ";
        //     sum += ntheta[i] * features[i];
        // }

        float sum = inner_product(ntheta.begin(), ntheta.end(), features.begin(), 0);
        return 1 - target * sum;
    }

    vector<float> sum(vector<float> &a, vector<float> &b)
    {
        vector<float> sum(nFeatures, 0.0);
        for (unsigned i = 0; i < nFeatures; i++)
        {
            // cout << ntheta[i] << "*" << features[i] << " ";
            sum[i] = a[i] + b[i];
        }
        return sum;
    }

    float _D_z(int i, float z)
    {
        vector<float> w = theta;
        w[i] = w[i] + z;
        float sum = 0.0;
        for (unsigned j = 0; j < mExamples; j++)
        {
            TrainingExample s = ts[j];
            float b = _b_func(w, s.getFeatures(), s.getTarget());
            if (b > 0.0)
            {
                sum = sum + b;
            }
        }
        return 0.5 * pow(w[i], 2) + C * pow(sum, 2);
    }

    bool _converge_condition(float a, float b, float z, float sigma = 0.01)
    {
        if ((a - b) <= (z * pow(sigma, 2)))
        {
            return true;
        }
        return false;
    }

    float _get_z(float d0, int i, vector<float> &ntheta, float beta = 0.5)
    {
        int k = 0;
        float lam = 1.0;
        float z_candidate = lam * d0;
        float D1 = _D_z(i, z_candidate);
        float D0 = _D_z(i, 0);
        while (!_converge_condition(D1, D0, z_candidate))
        {
            k = k + 1;
            lam = pow(beta, k);
            z_candidate = d0 * lam;
            D1 = _D_z(i, z_candidate);
            // cout << "loop" << k << endl << endl;
        }
        return z_candidate;
    }

    float _newton_direction(int i, float z, vector<float> &ntheta, vector<TrainingExample> &samples)
    {
        float sum1 = 0.0;
        float sum2 = 0.0;
        vector<float> w = ntheta;
        w[i] = w[i] + z;
        for (unsigned j = 0; j < mExamples; j++)
        {
            TrainingExample s = samples[j];
            float b = _b_func(w, s.getFeatures(), s.getTarget());
            if (b > 0.0)
            {
                float f = samples[j].getFeature(i);
                sum1 = sum1 + b * s.getTarget() * f;
                sum2 = sum2 + pow(f, 2);
            }
        }
        sum1 = ntheta[i] + z - 2 * C * sum1;
        sum2 = 1 + 2 * C * sum2;
        return sum1 / sum2;
    }

public:
    SVM(vector<TrainingExample> &examples)
    {
        nFeatures = 0;
        C = 1;
        ts = examples;
        mExamples = examples.size();
        if (mExamples > 0)
            nFeatures = examples[0].getFeatures().size();
        theta = vector<float>(nFeatures, 0.5);
        z = vector<float>(nFeatures, 0.0);
    }

    template <class RandomIt>
    void random_shuffle(RandomIt first, RandomIt last)
    {
        typedef typename std::iterator_traits<RandomIt>::difference_type diff_t;

        for (diff_t i = last - first - 1; i > 0; --i)
        {
            using std::swap;
            swap(first[i], first[std::rand() % (i + 1)]);
            // rand() % (i + 1) is not actually correct, because the generated number
            // is not uniformly distributed for most values of i. A correct implementation
            // will need to essentially reimplement C++11 std::uniform_int_distribution,
            // which is beyond the scope of this example.
        }
    }

    vector<float> coodinateDescent(int K = 100)
    {
        const float alpha = 0.0000001;
        const float eps = 0.00001;
        bool converge = false;
        int debug = 0;
        float diff = 0;
        diff = J();

        cout << "J(theta) = " << diff << endl
             << endl;

        for (unsigned t = 0; t < K; t++)
        {
            // cout << "Using example" << i << endl;
            vector<float> indices;

            if (t == 0)
            {
                for (unsigned i = 0; i < nFeatures; i++)
                {
                    // cout << ntheta[i] << "*" << features[i] << " ";
                    float d0 = 0.0 - _newton_direction(i, z[i], theta, ts);
                    float z_check = _get_z(d0, i, theta);
                    theta[i] = theta[i] + z_check;
                    indices.push_back(i);
                    cout << i << endl;
                    // if ((diff - J()) / abs(diff) >= 0.01)
                    // {
                    //     cout << "J(theta) = " << J() << endl << endl;
                    //     return theta;
                    // }
                }
            }
            else
            {
                random_shuffle(indices.begin(), indices.end());
                for (unsigned i = 0; i < nFeatures; i++)
                {
                    int pi = indices[i];
                    float d0 = 0.0 - _newton_direction(i, z[pi], theta, ts);
                    float z_check = _get_z(d0, i, theta);
                    theta[i] = theta[i] + z_check;
                    cout << i << endl;
                }
            }
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
