// Distributed under GNU General Public License (see license.txt for details).
//
//  Copyright (c) 2007 Shai Shalev-Shwartz.
//  All Rights Reserved.
//=============================================================================
// File Name: WeightVector.cc
// implements the methods of the WeightVector class
//=============================================================================
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sstream>

#include "simple_sparse_vec_hash.h"


class WeightVector
{

private:
  WeightVector(const WeightVector &); // disallowed

public:
  // Construct a Weight vector with dimension d
  WeightVector(uint dd) : d(dd), my_a(1.0), my_snorm(0.0), my_v(NULL)
  {
    my_v = new double[d];
    for (uint i = 0; i < d; ++i)
      my_v[i] = 0.0;
  }

  // Construct a Weight vector with dimension d from a file
  WeightVector(uint dd, std::ifstream &model_file) : d(dd), my_a(1.0), my_snorm(0.0), my_v(NULL)
  {
    my_v = new double[d];
    for (uint i = 0; i < d; ++i)
      my_v[i] = 0.0;
    unsigned int ind = 0;
    model_file >> ind;
    while (model_file.good())
    {
      char c;
      model_file >> c;
      double val;
      model_file >> val;
      my_v[ind] = val;
      my_snorm += val * val;
      model_file >> ind;
    }
  }

  // destructor
  ~WeightVector() { delete[] my_v; }

  // this *= s
  void scale(double s)
  {
    my_snorm *= (s * s);
    if (s != 0.0)
    {
      my_a *= s;
    }
    else
    {
      my_a = 1.0;
      for (uint i = 0; i < d; ++i)
        my_v[i] = 0.0;
    }
  }

  // this += s*x
  void add(simple_sparse_vector &x, double s);

  double getValue(int i);

  void addVec(float s, simple_sparse_vector &x);

  // this += s*x
  void add(WeightVector &x, double s);

  // this += s*x
  void add(float s);

  double dot(simple_sparse_vector &x);
  // this += s*x
  void add_dim(int i, double s);
  void assign(int i, double s);

  void print(std::ostream &os);

  double operator[](uint i)
  {
    if (i < d)
      return (my_v[i] * my_a);
    return 0.0;
  }

  uint dimension() { return d; }

  // ||this||^2
  double snorm()
  {
    return my_snorm;
  }

  // copy other
  void operator=(const WeightVector &other)
  {
    if (d != other.d)
    {
      std::cerr << "Assigning WeightVector of size "
                << other.d << " to a WeightVector of size "
                << d << " is not permitted" << std::endl;
      exit(EXIT_FAILURE);
    }
    my_a = other.my_a;
    my_snorm = other.my_snorm;
    for (uint i = 0; i < d; ++i)
      my_v[i] = other.my_v[i];
  }

  // make_my_a_one
  // use it for forcing my_a to be 1.0
  void make_my_a_one()
  {
    for (uint i = 0; i < d; ++i)
    {
      my_v[i] *= my_a;
    }
    my_a = 1.0;
  }

private:
  // The internal representation of w is as w = a*v where:
  uint d;
  double my_a;
  double my_snorm;
  double *my_v;
};

void WeightVector::add(simple_sparse_vector &x, double s)
{
  double pred = 0.0, norm_x = 0.0;
  for (simple_sparse_vector_iterator it = x.my_vec.begin();
       it != x.my_vec.end(); it++)
  {
    double val = (*it).second * s;
    norm_x += val * val;
    pred += 2.0 * my_v[(*it).first] * val;
    my_v[(*it).first] += (val / my_a);
  }
  my_snorm += norm_x + my_a * pred;
}

// this += s*x
void WeightVector::add(WeightVector &x, double s)
{
  my_snorm = 0.0;
  for (uint i = 0; i < d; ++i)
  {
    my_v[i] *= my_a;
    my_v[i] += (x[i] * s);
    my_snorm += my_v[i] * my_v[i];
  }
  my_a = 1.0;
}

// this += s*x
void WeightVector::add(float s)
{
  my_snorm = 0.0;
  for (uint i = 0; i < d; ++i)
  {
    my_v[i] += s;
    my_snorm += my_v[i] * my_v[i];
  }
  my_a = 1.0;
}

// this += s*x
void WeightVector::addVec(float s, simple_sparse_vector &x)
{
  for (simple_sparse_vector_iterator it = x.my_vec.begin();
       it != x.my_vec.end(); it++)
  {
    int index = (*it).first;
    my_v[index] += s * (*it).second;
  }
}

// this += s*x
double WeightVector::dot(simple_sparse_vector &x)
{
  double dot = 0.0;
  for (simple_sparse_vector_iterator it = x.my_vec.begin();
       it != x.my_vec.end(); it++)
  {
    int index = (*it).first;
    dot += my_v[index] * (*it).second;
  }
  return dot;
}

void WeightVector::add_dim(int i, double s)
{
  // my_snorm = 0.0;
  my_v[i] += s;
}

double WeightVector::getValue(int i)
{
  // my_snorm = 0.0;
  return my_v[i];
}

void WeightVector::assign(int i, double s)
{
  // my_snorm = 0.0;
  my_v[i] = s;
}

void WeightVector::print(std::ostream &os)
{

  for (uint i = 0; i < d; ++i)
  {
    if (my_v[i] != 0.0)
      os << i << ":" << (my_v[i] * my_a) << " ";
  }
  os << std::endl;
}

//--------------------------------------------------------------------------
double operator*(simple_sparse_vector &u, WeightVector &v)
{
  double outcome = 0.0;
  for (simple_sparse_vector_iterator it = u.my_vec.begin();
       it != u.my_vec.end(); it++)
  {
    outcome += ((*it).second * v[(*it).first]);
  }
  return outcome;
}

//-----------------------------------------------------------------------------
double operator*(WeightVector &v, simple_sparse_vector &u)
{
  return (u * v);
}
