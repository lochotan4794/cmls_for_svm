// Distributed under GNU General Public License (see license.txt for details).
//
//  Copyright (c) 2007 Shai Shalev-Shwartz.
//  All Rights Reserved.
//
// Main for the PEGASOS algorithm.
//
//
//
// Pegasos: Fly high in the primal
// main file
//
//

//*****************************************************************************
// Included Files
//*****************************************************************************
#include "simple_sparse_vec_hash.h"
#include "greedy_cd.h"
// #include "WeightVector.h"

int main(int argc, char **argv)
{

  // -------------------------------------------------------------
  // ---------------------- Parse Command Line -------------------
  // -------------------------------------------------------------

  std::string data_filename = "B_train.dat";
  // std::string model_filename = "noModelFile";
  
  std::string model_filename = "model";
  double lambda = 1.0;
  int max_iter = 10;
  int exam_per_iter = 1;
  // -------------------------------------------------------------
  // ---------------------- read the data ------------------------
  // -------------------------------------------------------------
  uint dimension = 0;
  std::vector<simple_sparse_vector> Dataset;
  std::vector<int> Labels;
  long readingTime;
  Lasso model = Lasso();
  model.ReadData(data_filename, Dataset, Labels, dimension, readingTime);

  uint testDimension = 0;
  std::vector<simple_sparse_vector> testDataset;
  std::vector<int> testLabels;
  long testReadingTime;
  std::string test_filename = "B_test.dat";
  model.ReadData(test_filename, testDataset, testLabels, testDimension, testReadingTime);

  std::cerr << readingTime + testReadingTime << " = Time for reading the data" << std::endl;

  // choose a random seed
  srand(time(NULL));

  // -------------------------------------------------------------
  // ---------------------- Experiments mode ---------------------
  // -------------------------------------------------------------
  // if (experiments_file != "noExperimentsFile")
  // {
  //   run_experiments(experiments_file, Dataset, Labels, dimension,
  //                   testDataset, testLabels);
  //   return (EXIT_SUCCESS);
  // }

  // -------------------------------------------------------------
  // ---------------------- Main Learning function ---------------
  // -------------------------------------------------------------
  long trainTime, calc_obj_time;
  double obj_value, norm_value, loss_value, zero_one_error, test_loss, test_error;
  /*
  Learn(Dataset,Labels,dimension,testDataset,testLabels,
  lambda,max_iter,exam_per_iter,
  num_iter_to_avg,model_filename,
  trainTime,calc_obj_time,obj_value,norm_value,loss_value,zero_one_error,
  test_loss,test_error,
  0,0.0,0,0.0);
  uint num_example_to_validate = max_iter/10;
  LearnAndValidate(Dataset,Labels,dimension,testDataset,testLabels,
       lambda,max_iter,exam_per_iter,
       num_example_to_validate,model_filename,
       trainTime,calc_obj_time,obj_value,norm_value,
       loss_value,zero_one_error,
       test_loss,test_error,
       0,0.0,0,0.0);
  */
  double eps = 0.1;
  model.LearnReturnLast(Dataset, Labels, dimension, testDataset, testLabels,
                        lambda, max_iter, eps, exam_per_iter,
                        model_filename,
                        trainTime, calc_obj_time, obj_value, norm_value,
                        loss_value, zero_one_error,
                        test_loss, test_error,
                        0, 0.0, 0, 0.0);

  // -------------------------------------------------------------
  // ---------------------- Print Results ------------------------
  // -------------------------------------------------------------
  std::cout << trainTime << " = Time for training\n"
            << calc_obj_time << " = Time for calculate objective\n"
            << norm_value << " = Norm of solution\n"
            << loss_value << " = avg Loss of solution\n"
            << zero_one_error << " = avg zero-one error of solution\n"
            << obj_value << " = primal objective of solution\n"
            << test_loss << " = avg Loss over test\n"
            << test_error << " = avg zero-one error over test\n"
            << std::endl;

  return (EXIT_SUCCESS);
}



