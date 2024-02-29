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
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <iomanip>

#include "simple_sparse_vec_hash.h"
#include "WeightVector.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::shuffle;
using std::swap;
using std::vector;

class SVM
{
private:
    float C;
    int nExamples;

    // help function for getting runtime
    long get_runtime(void)
    {
        clock_t start;
        start = clock();
        return ((long)((double)start / (double)CLOCKS_PER_SEC));
    }

    float kernel(simple_sparse_vector &p1, simple_sparse_vector &p2)
    {
        return dotProduct(p1, p2);
    }

    float dotProduct(simple_sparse_vector &v1, simple_sparse_vector &v2)
    {

        int i = 0, j = 0;
        float result = 0.0;

        while (i < v1.my_vec.size() && j < v2.my_vec.size())
        {
            if (v1.my_vec[i].first < v2.my_vec[j].first)
            {
                i++;
            }
            else if (v1.my_vec[i].first > v2.my_vec[j].first)
            {
                j++;
            }
            else
            {
                result += (v1.my_vec[i].second * v2.my_vec[j].second);
                i++;
                j++;
            }
        }

        return result;
    }

    float obj(int i1, int i2, float a1, float a2, vector<simple_sparse_vector> &Dataset, vector<int> &Labels, WeightVector &W, vector<float> &alpha, float **kernelMatrix)
    {
        float return_value = 0;

        int y1 = Labels[i1];
        int y2 = Labels[i2];
        simple_sparse_vector s1 = Dataset[i1];
        simple_sparse_vector s2 = Dataset[i2];
        float k11 = kernelMatrix[i1][i1];
        float k12 = kernelMatrix[i1][i2];
        float k22 = kernelMatrix[i2][i2];
        float k21 = kernelMatrix[i2][i1];

        return_value = k11 * alpha[i1] * a1 * y1 * y1 + k12 * alpha[i1] * a2 * y1 * y2 + k21 * a2 * alpha[i1] * y2 * y1 + k22 * a2 * a2 * y2 * y2;

        float sum = alpha[i1] + a2;

        return sum - 0.5 * return_value;
    }

    int heuristic(int i2, vector<float> &E, int mExamples)
    {
        float maxDiff = 0;
        int index = 0;
        if (E[i2] > 0)
        {
            for (unsigned j = 0; j < mExamples; j++)
            {
                if (E[j] < maxDiff && i2 != j)
                {
                    index = j;
                    maxDiff = E[j];
                }
            }
        }
        else
        {
            for (unsigned j = 0; j < mExamples; j++)
            {
                if (E[j] > maxDiff && i2 != j)
                {
                    index = j;
                    maxDiff = E[j];
                }
            }
        }
        return index;
    }


    float takeStep(int i1, int i2, float eps, float mu, vector<simple_sparse_vector> &Dataset, vector<int> &Labels, WeightVector &W, vector<float> &E, vector<float> &alpha, float &b, float &C, int &num_ilterations, float **kernelMatrix)
    {
        simple_sparse_vector s1 = Dataset[i1];
        simple_sparse_vector s2 = Dataset[i2];
        float y1 = Labels[i1];
        float y2 = Labels[i2];
        float alph1 = alpha[i1];
        float alph2 = alpha[i2];
        float E1 = out(s1, W, b) - y1;
        // float E2 = out(s2, W, b) - y2;
        float s = y1 * y2;
        float L;
        float H;
        float a2 = 0;
        float a1 = 0;
        float b1 = 0;
        float b2 = 0;
        float Lobj;
        float Hobj;
        float b_n;
        if (y1 != y2)
        {
            L = 0.0 > (alph2 - alph1) ? 0.0 : (alph2 - alph1);
            H = min(C, C + alph2 - alph1);
        }
        else
        {
            L = 0.0 > (alph2 + alph1 - C) ? 0.0 : (alph2 + alph1 - C);
            H = min(C, alph2 + alph1);
        }
        if (i1 == i2)
        {
            return 0;
        }
        if (L == H)
        {
            return 0;
        }
        float k11 = kernelMatrix[i1][i1];
        float k12 = kernelMatrix[i1][i2];
        float k22 = kernelMatrix[i2][i2];
        float eta = -k11 - k22 + 2 * k12;
        if (eta < 0)
        {
            a2 = alph2 - y2 * (E1 - E[i2]) / eta;
            if (a2 < L)
                a2 = L;
            else if (a2 > H)
                a2 = H;
        }
        else
        {
            Lobj = obj(i1, i2, 0, L, Dataset, Labels, W, alpha, kernelMatrix);
            Hobj = obj(i1, i2, 0, H, Dataset, Labels, W, alpha, kernelMatrix);
            if (Lobj < Hobj - eps)
            {
                a2 = H;
            }
            else if (Lobj > Hobj + eps)
            {
                a2 = L;
            }
            else
            {
                a2 = alph2;
            }
        }
        if (a2 < 1e-8)
        {
            a2 = 0;
        }
        else if (a2 > C - 1e-8)
        {
            a2 = C;
        }
        if (abs(a2 - alph2) < eps * (a2 + alph2 + eps))
        {
            return 0;
        }
        a1 = alph1 + s * (alph2 - a2);
        float b1_n = E1 + y1 * (a1 - alph1) * kernel(s1, s2) + y2 * (a2 - alph2) * kernel(s1, s2) + b;
        float b2_n = E[i2] + y1 * (a1 - alph1) * kernel(s1, s2) + y2 * (a2 - alph2) * kernel(s1, s2) + b;
        bool bound = false;
        if (a1 != 0 && a1 != C)
        {
            b_n = b1;
        }
        else if (a2 != 0 && a2 != C)
        {
            b_n = b2;
        }
        else if ((a1 == 0 || a1 == C) && (a2 == 0 || a2 == C) && (L != H))
        {
            b_n = 0.5 * (b1_n + b2_n);
            bound = true;
        }

        W.addVec(y1 * (a1 - alph1), s1);
        W.addVec(y2 * (a2 - alph2), s2);
        static int num = 0;

        num = num + 1;
        cout << num << endl;

        for (unsigned i = 0; i < Labels.size(); i++)
        {
            if (i != i1 && i != i2)
                E[i] = E[i] + y1 * (a1 - alph1) * kernelMatrix[i1][i] + y2 * (a2 - alph2) * kernelMatrix[i2][i] + b - b_n;
        }

        b = b_n;

        alpha[i1] = a1;
        alpha[i2] = a2;

        return 1;
    }

    float examineExample(int i2, vector<simple_sparse_vector> &Dataset, vector<int> &Labels, WeightVector &W, vector<float> &E, vector<float> &alpha, float tol, float &b, float &C, int &num_ilterations, float **kernelMatrix)
    {
        simple_sparse_vector s = Dataset[i2];
        float y2 = Labels[i2];
        float alph2 = alpha[i2];
        float E2 = out(s, W, b) - y2;
        E[i2] = E2;
        float r2 = E2 * y2;
        int mExamples = Dataset.size();

        if ((r2 < -tol && alph2 < C) || (r2 > tol && alph2 > 0))
        {
            int nn = 0;
            vector<int> Index(mExamples, 0);
            for (unsigned j = 0; j < mExamples; j++)
            {
                if (alpha[j] != 0 && alpha[j] != C)
                {
                    nn++;
                    Index[j] = 1;
                }
                else
                {
                    Index[j] = 0;
                }
            }
            if (nn > 1)
            {
                int i1 = heuristic(i2, E, mExamples);
                if (takeStep(i1, i2, 0.001, 0.01, Dataset, Labels, W, E, alpha, b, C, num_ilterations, kernelMatrix))
                    return 1;
            }

            std::random_device rd;
            std::mt19937 g(rd());

            shuffle(Index.begin(), Index.end(), g);

            for (auto j = Index.begin(); j != Index.end(); ++j)
            {

                int i1 = *j;
                if (takeStep(i1, i2, 0.001, 0.01, Dataset, Labels, W, E, alpha, b, C, num_ilterations, kernelMatrix))
                    return 1;
            }

            std::srand(unsigned(std::time(0)));
            std::vector<int> myvector;

            // set some values:
            for (int i = 1; i < mExamples; ++i)
                myvector.push_back(i); // 1 2 3 4 5 6 7 8 9

            // using built-in random generator:
            shuffle(myvector.begin(), myvector.end(), g);

            for (auto j = myvector.begin(); j != myvector.end(); ++j)
            {
                if (*j != i2)
                {
                    if (takeStep(*j, i2, 0.001, 0.01, Dataset, Labels, W, E, alpha, b, C, num_ilterations, kernelMatrix))
                        return 1;
                }
            }
        }
        return 0;
    }

    float dot(vector<float> &ntheta, vector<float> &features, float C = 1)
    {
        float sum = 0.0;
        for (unsigned i = 0; i < features.size(); i++)
        {
            sum += ntheta[i] * features[i];
        }
        return sum;
    }

    float out(simple_sparse_vector &features, WeightVector &W, float &b)
    {
        float sum = 0.0;
        sum = W.dot(features) - b;
        return sum;
    }

public:
    SVM()
    {
        C = 10;
    }

    // ------------------------------------------------------------//
    // ---------------- READING DATA ------------------------------//
    // ------------------------------------------------------------//
    void ReadData( // input
        std::string &data_filename,
        // output
        std::vector<simple_sparse_vector> &Dataset,
        std::vector<int> &Labels,
        uint &dimension,
        long &readingTime)
    {

        dimension = 0;

        // Start a timer
        long startTime = get_runtime();

        // OPEN DATA FILE
        // =========================
        std::ifstream data_file(data_filename.c_str());
        if (!data_file.good())
        {
            std::cerr << "error w/ " << data_filename << std::endl;
            exit(EXIT_FAILURE);
        }

        // Read SVM-Light data file
        // ========================
        int num_examples = 0;
        std::string buf;
        while (getline(data_file, buf))
        {
            // ignore lines which begin with #
            if (buf[0] == '#')
                continue;
            // Erase what comes after #
            size_t pos = buf.find('#');
            if (pos < buf.size())
            {
                buf.erase(pos);
            }
            // replace ':' with white space
            int n = 0;
            for (size_t pos = 0; pos < buf.size(); ++pos)
                if (buf[pos] == ':')
                {
                    n++;
                    buf[pos] = ' ';
                }
            // read from the string
            std::istringstream is(buf);
            int label = 0;
            is >> label;
            // if (label != 1 && label != -1)
            // {
            //     std::cerr << "Error reading SVM-light format. Abort." << std::endl;
            //     exit(EXIT_FAILURE);
            // }
            if (label == 0)
            {
                label = -1;
            }

            Labels.push_back(label);
            simple_sparse_vector instance(is, n);
            Dataset.push_back(instance);
            num_examples++;
            uint cur_max_ind = instance.max_index() + 1;
            if (cur_max_ind > dimension)
                dimension = cur_max_ind;
        }

        data_file.close();

#ifdef nodef
        std::cerr << "num_examples = " << num_examples
                  << " dimension = " << dimension
                  << " Dataset.size = " << Dataset.size()
                  << " Labels.size = " << Labels.size() << std::endl;
#endif

        // update timeline
        readingTime = get_runtime() - startTime;
    }

    float **getKernel(int N, vector<simple_sparse_vector> &Dataset)
    {
        float **arr = new float *[N];
        for (int i = 0; i < N; ++i)
        {
            arr[i] = new float[N];
            for (int j = 0; j < N; ++j)
            {
                simple_sparse_vector s1 = Dataset[i];
                simple_sparse_vector s2 = Dataset[j];

                arr[i][j] = kernel(s1, s2);
            }
        }
        return arr;
    }

    void LearnReturnLast( // Input variables
        vector<simple_sparse_vector> &Dataset,
        vector<int> &Labels,
        uint dimension,
        std::vector<simple_sparse_vector> &testDataset,
        std::vector<int> &testLabels,
        double lambda, int max_iter, double eps,
        int exam_per_iter,
        std::string &model_filename,
        // Output variables
        long &train_time, long &calc_obj_time,
        double &obj_value, double &norm_value,
        double &loss_value, double &zero_one_error,
        double &test_loss, double &test_error,
        // additional parameters
        int eta_rule_type, double eta_constant,
        int projection_rule, double projection_constant, int num_ilterations)
    {

        long startTime = get_runtime();
        long endTime;
        uint num_examples = Labels.size();
        nExamples = num_examples;

        // Initialization of classification vector
        WeightVector W(dimension);
        vector<float> E(num_examples, 0.0);

        // ---------------- Main Loop -------------------
        cout << "dimension: " << num_examples << endl;

        float numChanged = 0;
        float examineAll = 1;
        vector<float> alpha(num_examples, 0.0);
        float b = 0.0;
        int count = 0;
        W.print(std::cout);
        float C = 1;
        startTime = get_runtime();

        float **kernelMatrix = getKernel(nExamples, Dataset);
        std::clock_t c_start = std::clock();

        while ((numChanged > 0 || examineAll) && count < 1)
        {
            numChanged = 0;
            if (examineAll)
            {
                for (unsigned i = 0; i < num_examples; i++)
                {
                    numChanged += examineExample(i, Dataset, Labels, W, E, alpha, 0.02, b, C, num_ilterations, kernelMatrix);
                }
            }
            else
            {
                for (unsigned i = 0; i < num_examples; i++)
                {
                    if (alpha[i] != 0.0 && alpha[i] != C)
                    {
                        cout << i << endl;
                        simple_sparse_vector s1 = Dataset[i];
                        numChanged += examineExample(i, Dataset, Labels, W, E, alpha, 0.02, b, C, num_ilterations, kernelMatrix);
                    }
                }
            }
            if (examineAll == 1)
            {
                examineAll = 0;
            }
            else if (numChanged == 0)
            {
                examineAll = 1;
            }
            count += 1;
            cout << count << endl;
        }
        // update timeline
        std::clock_t c_end = std::clock();
        endTime = get_runtime();
        train_time = endTime - startTime;
        cout << train_time << endl;
        float time_elapsed_ms = 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC;
        W.print(std::cout);
        std::cout << "CPU time used: " << time_elapsed_ms << " ms\n";
        std::cout << "num_ilterations: " << num_ilterations << " ms\n";

        // Calculate objective value
        norm_value = W.snorm();
        obj_value = norm_value * lambda / 2.0;
        loss_value = 0.0;
        zero_one_error = 0.0;
        for (uint i = 0; i < Dataset.size(); ++i)
        {
            double cur_loss = 1 - Labels[i] * (W * Dataset[i] - b);
            if (cur_loss < 0.0)
                cur_loss = 0.0;
            loss_value += cur_loss / num_examples;
            obj_value += cur_loss / num_examples;
            if (cur_loss >= 1.0)
                zero_one_error += 1.0 / num_examples;
        }

        endTime = get_runtime();
        calc_obj_time = endTime - startTime;

        // Calculate test_loss and test_error
        test_loss = 0.0;
        test_error = 0.0;
        for (uint i = 0; i < testDataset.size(); ++i)
        {
            double cur_loss = 1 - testLabels[i] * (W * testDataset[i] - b);
            if (cur_loss < 0.0)
                cur_loss = 0.0;
            test_loss += cur_loss;
            if (cur_loss >= 1.0)
                test_error += 1.0;
        }
        if (testDataset.size() != 0)
        {
            test_loss /= testDataset.size();
            test_error /= testDataset.size();
        }

        // finally, print the model to the model_file
        if (model_filename != "noModelFile")
        {
            std::ofstream model_file(model_filename.c_str());
            if (!model_file.good())
            {
                std::cerr << "error w/ " << model_filename << std::endl;
                exit(EXIT_FAILURE);
            }
            W.print(model_file);
            model_file.close();
        }
    }
};

#endif
