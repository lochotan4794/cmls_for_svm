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
#include "simple_sparse_vec_hash.h"
#include "WeightVector.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::swap;

class SVM
{
private:
    // WeightVector &W;
    // std::vector<simple_sparse_vector> z;
    // std::vector<simple_sparse_vector> &Dataset;
    // std::vector<int> &Labels;
    // unsigned int mExamples, nFeatures;
    double C;
    int K;

    // help function for getting runtime
    long get_runtime(void)
    {
        clock_t start;
        start = clock();
        return ((long)((double)start / (double)CLOCKS_PER_SEC));
    }

    double D_z(int C, int i, float z, std::vector<simple_sparse_vector> &Dataset, std::vector<int> &Label, WeightVector &W)
    {
        WeightVector &theta = W;
        theta.add_dim(i, z);
        float sum = 0.0;
        for (unsigned t = 0; t < Dataset.size(); t++)
        {
            simple_sparse_vector s = Dataset[t];
            float b = b_func(W, s, Label[t]);
            if (b > 0.0)
            {
                sum = sum + b;
            }
        }
        return 0.5 * pow(theta.getValue(i), 2) + C * pow(sum, 2);
    }

    bool _converge_condition(float a, float b, float z, float sigma = 0.01)
    {
        if ((a - b) <= (z * pow(sigma, 2)))
        {
            return true;
        }
        return false;
    }

    double get_z(int C, float d0, int i, WeightVector &W, double beta, std::vector<simple_sparse_vector> &Dataset, std::vector<int> &Label)
    {
        int k = 0;
        float lam = 1.0;
        float z_candidate = lam * d0;
        float D1 = D_z(C, i, z_candidate, Dataset, Label, W);
        float D0 = D_z(C, i, 0, Dataset, Label, W);

        while (!_converge_condition(D1, D0, z_candidate) && z_candidate > 0.0)
        {
            k = k + 1;
            lam = pow(beta, k);
            z_candidate = d0 * lam;
            cout << z_candidate << endl;
            D1 = D_z(C, i, z_candidate, Dataset, Label, W);
            cout << "loop" << k << endl
                 << endl;
        }
        return z_candidate;
    }

    double b_func(WeightVector &W, simple_sparse_vector &features, int label)
    {
        double dot_value = W.dot(features);
        return 1 - label * dot_value;
    }

    float newton_direction(int C, int i, float z, WeightVector &W, std::vector<simple_sparse_vector> &Dataset, std::vector<int> &Label)
    {
        float sum1 = 0.0;
        float sum2 = 0.0;
        WeightVector &theta = W;
        theta.add_dim(i, z);
        for (unsigned j = 0; j < Dataset.size(); j++)
        {
            simple_sparse_vector s = Dataset[j];
            float b = b_func(theta, s, Label[j]);
            if (b > 0.0)
            {
                float f = s.get_second(i);
                sum1 = sum1 + b * Label[j] * f;
                sum2 = sum2 + pow(f, 2);
            }
        }
        sum1 = theta.getValue(i) + z - 2 * C * sum1;
        sum2 = 1 + 2 * C * sum2;
        return sum1 / sum2;
    }

    double CK(double r, double c)
    {
        if (r <= 1)
        {
            return 2 * (r - 1);
        }
        else
        {
            return 2 * c * (r - 1);
        }
    }

    double FK(double x, double y, double c)
    {
        if (x <= 1 + abs(y))
        {
            return 2;
        }
        else
        {
            return 2 * c;
        }
    }

public:
    SVM(double C = 1.0, int K = 100)
    {
        C = C;
        K = K;
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
            if (label != 1 && label != -1)
            {
                std::cerr << "Error reading SVM-light format. Abort." << std::endl;
                exit(EXIT_FAILURE);
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

    void LearnReturnLast( // Input variables
        std::vector<simple_sparse_vector> &Dataset,
        std::vector<int> &Labels,
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
        int projection_rule, double projection_constant)
    {

        long startTime = get_runtime();
        long endTime;
        int K = 100;
        double z = 0;
        uint num_examples = Labels.size();
        std::vector<double> delta(dimension, 10);

        // Initialization of classification vector
        WeightVector W(dimension);
        std::vector<double> R(num_examples, 0.0);

        // ---------------- Main Loop -------------------
        cout << "dimension: " << dimension << endl;

        for (int t = 1; t < K; ++t)
        {
            double c = max(0, 1 - t /50);

            for (int j = 0; j < dimension; ++j)
            {
                double k1, k2 = 0;
                for (int i = 0; i < num_examples; ++i)
                {
                    k1 += Dataset[i].get_second(j) * Labels[i] * CK(R.at(i), c) + 2 * lambda * num_examples * W.getValue(j);
                    k2 += FK(R.at(i), delta.at(j) * Dataset[i].get_second(j), c) * pow(Dataset[i].get_second(j), 2) + 2 * lambda * num_examples;
                }
                double deltaV = -k1 / k2;
                double deltaW = min(max(deltaV, -delta[j]), delta[j]);
                for (int i = 0; i < num_examples; ++i)
                {
                    R[i] = R[i] + deltaW * Dataset[i].get_second(j) * Labels[i];
                }
                W.add_dim(j, deltaW);
                delta[j] = 2 * abs(deltaW) + eps;

                // if (j % 1000 == 0)
                //     cout << "j= " << j << endl;
            }
            cout << "t= " << t << endl;
        }

        // update timeline
        endTime = get_runtime();
        train_time = endTime - startTime;
        startTime = get_runtime();

        // Calculate objective value
        norm_value = W.snorm();
        obj_value = norm_value * lambda / 2.0;
        loss_value = 0.0;
        zero_one_error = 0.0;
        for (uint i = 0; i < Dataset.size(); ++i)
        {
            double cur_loss = 1 - Labels[i] * (W * Dataset[i]);
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
            double cur_loss = 1 - testLabels[i] * (W * testDataset[i]);
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
