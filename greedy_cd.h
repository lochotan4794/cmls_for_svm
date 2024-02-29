
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
#include <errno.h>

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::swap;

class Lasso
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

    double sign(double x)
    {
        if (x > 0)
        {
            return 1;
        }
        else
        {
            return -1;
        }
    }

    double shrinkage(double x, double param)
    {
        if (abs(x) > param)
        {
            return x - sign(x) * param;
        }
        else
        {
            return 0.0;
        }
    }

public:
    Lasso(double C = 1.0, int K = 100)
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
            if (label != 1 && label != 0)
            {
                cout << label << endl;
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
        int K = 10;
        double z = 0;
        double L = 1.0;
        uint num_examples = Labels.size();
        std::vector<double> delta(dimension, 10);

        // Initialization of classification vector
        WeightVector W(dimension);
        std::vector<double> R(num_examples, 0.0);

        // ---------------- Main Loop -------------------
        cout << "dimension: " << dimension << endl;

        for (int t = 0; t < K; ++t)
        {
            int index = 0;
            double Ji;
            double tmp;
            for (uint i = 0; i < dimension; ++i)
            {
                double cur_loss = 0.0;
                for (uint j = 0; j < Dataset.size(); ++j)
                {
                    double w;
                    w = W.getValue(i);
                    cur_loss += w * (Labels[j] - w * Dataset[j].get_second(i));
                }
                cur_loss = cur_loss / Dataset.size();
                if (W[i] == 0)
                {
                    tmp = shrinkage(cur_loss, lambda);
                }
                else
                {
                    tmp = cur_loss + sign(W[i]) * lambda;
                }
                if (tmp > Ji)
                {
                    index = i;
                    Ji = tmp;
                }
            }
            double alpha_plus = shrinkage(W.getValue(index) - (1 / L) * Ji, lambda / L);
            cout << alpha_plus << endl;
            if ((alpha_plus * W[index]) >= 0)
            {
                cout << "update" << endl;
                W.assign(index, alpha_plus);
            }
            else
            {
                W.assign(index, 0.0);
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
            double cur_loss = Labels[i] - (W * Dataset[i]);
            // if (cur_loss < 0.0)
            //     cur_loss = 0.0;
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
            double cur_loss = testLabels[i] - (W * testDataset[i]);
            // if (cur_loss < 0.0)
            //     cur_loss = 0.0;
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
