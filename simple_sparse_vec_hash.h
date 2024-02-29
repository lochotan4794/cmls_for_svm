// Distributed under GNU General Public License (see license.txt for details).
//
//  Copyright (c) 2007 Shai Shalev-Shwartz.
//  All Rights Reserved.
//==============================================================================
// File Name: simple_sparse_vec_hash.h
// Written by: Shai Shalev-Shwartz (28.01.07)
// implement a very simple hash table and sparse vector based on stl vector 
// and the modulu function
//==============================================================================

#ifndef _SHAI_SIMPLE_HASH_TABLE
#define _SHAI_SIMPLE_HASH_TABLE

//*****************************************************************************
// Included Files
//*****************************************************************************
#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <sstream>


class IndexValuePair {
 public:
  IndexValuePair() { }
  IndexValuePair(uint i, double v) : first(i), second(v) { }
  uint first;
  double second;
};
 
typedef std::vector<IndexValuePair>::iterator simple_sparse_vector_iterator;


/** Implements a simple sparse vector based on stl vector.
    @author Shai Shalev-Shwartz (shais@cs.huji.ac.il)
*/
class simple_sparse_vector {
 public:

  /** Default Constructor. Allocates an all zeros sparse vector
  */
  simple_sparse_vector() {  }

  /** Constructor. Read a line from a file of pairs of string and value.
      Construct a vector from this line.
      @param is The input stream to read from.
  */
  simple_sparse_vector(std::istream& is);

  /** Constructor. Construct a vector from an istringstream.
      @param is The input stringstream to read from.
      @param n The number of elements to read
  */
  simple_sparse_vector(std::istringstream& is,int n);
  

  // Multiply a sparse vector by a scalar s
  void scale(double s);
 
  /** Returns the squared l_2 norm of the vector
      @return the squared l_2 norm of the vector
  */
  double snorm();

  double dot(std::vector<double> x);

  /** Convert the vector to a binary vector that just indicates 
      which elements are non-zero
  */
  void make_binary();

  double get_second(int i);

  // return the maximal non-zero index
  uint max_index();

  // Print the content of a sparse instance to an output stream
  void print(std::ostream& os);

  // Zero all indexed elements of a sparse vector
  void zero();

  std::vector<IndexValuePair> my_vec;

};

/*---------------------------------------------------------------------------*/
simple_sparse_vector::simple_sparse_vector(std::istream& is) : my_vec() {
  
  // read the number of elements
  int n = 0;
  is >> n;

  // read the elements
  for (int i=0; i<n ; ++i) {

    // read the pair (key,val)
    uint key;
    is >> key;
    double val;
    is >> val;
    // insert to the map
    my_vec.push_back(IndexValuePair(key,val));

  }

}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
simple_sparse_vector::simple_sparse_vector(std::istringstream& is,int n) : my_vec() {
  
  // read the elements
  for (int i=0; i<n ; ++i) {

    // read the pair (key,val)
    uint key;
    is >> key;
    double val;
    is >> val;

    // insert to the map
    my_vec.push_back(IndexValuePair(key,val));

  }
}

/*---------------------------------------------------------------------------*/
double simple_sparse_vector::dot(std::vector<double> x){
  double res = 0.0;
  for( simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
        res += (*it).second * x[(*it).first];
  }
  return res;
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
void simple_sparse_vector::scale(double s) {

  for( simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
    (*it).second *= s;
  }
}

double simple_sparse_vector::get_second(int i) {

  for( simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
    if ((*it).first == i)
    {
      return (*it).second;
    }
  }
  return 0.0;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
uint simple_sparse_vector::max_index() {

  uint d=0;
  for( simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
    if ((*it).first > d) d = (*it).first;
  }
  return d;
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
double simple_sparse_vector::snorm() {

  double output = 0.0;
  for(simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
    double tmp = (*it).second;
    output += tmp*tmp;
  }

  return(output);
}
/*---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*/
void simple_sparse_vector::make_binary() {

  for(simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
    (*it).second = 1.0;
  }

}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void simple_sparse_vector::print(std::ostream& os) {

  for(simple_sparse_vector_iterator it = my_vec.begin(); 
      it != my_vec.end(); it++) {
    os << "(" << (*it).first << "," << (*it).second << ") "; 
  }
  os << std::endl;
}
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
void simple_sparse_vector::zero() {
  my_vec.clear();
}
/*---------------------------------------------------------------------------*/


#ifdef nodef


typedef std::vector< std::vector<IndexValuePair> >::iterator simple_hash_table_iterator;



/** Implements a simple hash table based on stl vector.
    The hash function is the modulu %
    @author Shai Shalev-Shwartz (shais@cs.huji.ac.il)
*/
class simple_hash_table {
 public:

  /** Default Constructor. Allocates an all zeros hash table
  */
  simple_hash_table() : my_vec(1087) {  }

  /** Constructor. Allocates an all zeros hash table with m entries
  */
  simple_hash_table(uint m) : my_vec(m) {  }


  /** retreive the value of the i'th element of the hash
      @param i the index to retrieve
  */
  double get(uint i);

  /** reference operator. Retrieves the i'th element of the hash and
      create it if it doesn't exist.
      @param i the index to retrieve
  */
  double& get_ref(uint i);

  // Multiply a sparse vector by a scalar s
  void scale(double s);

  // this = a * this + b * other;
  void scale_and_add(simple_sparse_vector& other, double a, double b);

  // this = a * this + b * other;
  void scale_and_add(simple_hash_table& other, double a, double b);

  // this += b*other;
  void add(simple_hash_table& other, double b);

  // this += b*other;
  void add(simple_sparse_vector& other, double b);


  /** Returns the squared l_2 norm of the vector
      @return the squared l_2 norm of the vector
  */
  double snorm();

  // Print the content of a sparse instance to an output stream
  void print(std::ostream& os);

  // Zero all indexed elements of a sparse vector
  void zero();


  std::vector< std::vector<IndexValuePair> > my_vec;

};


//-----------------------------------------------------------------------------
/** Operator * for vector-vector multiplication
    @param u A reference to a simple_sparse_vector
    @param v A reference to a simple_hash_table
    @return The product (double)
*/
double operator* (simple_sparse_vector& u, simple_hash_table& v);

//-----------------------------------------------------------------------------
/** Operator * for vector-vector multiplication
    @param u A reference to a simple_hash_table
    @param v A reference to a simple_sparse_vector
    @return The product (double)
*/
double operator* (simple_hash_table& u, simple_sparse_vector& v);

#endif

#endif
