#ifndef _VECTOR_H_
#define _VECTOR_H_
#include <concepts>
#include <vector>
#include <cblas.h>
#include <algorithm>
#include <iostream>

template<std::floating_point T>
class Vector
{
public:////////////////////////////////////////////////Constructors and Initialization
  Vector() : data_() {};
  explicit Vector(size_t size) : data_(size, 0.0){};/// Vector vec(3);
  Vector(const std::initializer_list<T>& init) : data_(init) {};/// Vector vec = [1.0, 2.0, 3.0];
  Vector(const Vector<T>& vec): data_(vec.data_){};/// copy constructor

  /// resize the Vector  
  void reinit(size_t size){
    ///If size is greater than current one, more elements will be
    ///added with the initial value 0, with old elements unchanged. If
    ///size is same to the current one, nothing would be done. If size
    ///is less than the current one, the tail will be just removed for
    ///fit the size.
    data_.resize(size);
  }

  void reinit(size_t size, bool fast){
    if (fast){
      data_.resize(size);
    }
    else {
      data_.resize(size);
      std::fill(data_.begin(), data_.end(), 0.);
    }
  }
  
  
public://////////////////////////////////////////////// Element Access and Modification
  // Element access
  T& operator()(size_t index) {
    return data_[index];
  }

  const T& operator()(size_t index) const {
    return data_[index];
  }

  /// Assignment operators
  Vector &operator = (const T val) {
    std::fill(data_.begin(), data_.end(), val);
    return *this;
  }  
  

public:///////////////////////////////////////////////// Arithmetic Operations
  // Scalar multiplication: this = scalar * this
  void scale(T scalar) {
    if constexpr (std::is_same_v<T, float>){
      cblas_sscal(data_.size(), scalar, data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_dscal(data_.size(), scalar, data_.data(), 1);
    }
  }
  
  void sadd (const T s, const Vector<T> &V){
    ///Scaling and simple vector addition,
    ///i.e. *<tt>*this = s*(*this)+V</tt>.
    
    if (data_.size() != V.data_.size()) {
      throw std::invalid_argument("Vectors must be of the same size");
    }

    // First, scale this vector by s
    if constexpr (std::is_same_v<T, float>){
      cblas_sscal(data_.size(), s, data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_dscal(data_.size(), s, data_.data(), 1);
    }

    // Then, add V to this vector
    if constexpr (std::is_same_v<T, float>){
      cblas_saxpy(data_.size(), 1.0, V.data_.data(), 1, data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_daxpy(data_.size(), 1.0, V.data_.data(), 1, data_.data(), 1);
    }
 
  }

  void add(T alpha, const Vector<T>& v) {
    /// Vector addition: this = this + v
    
    if (data_.size() != v.data_.size()) {
      throw std::invalid_argument("Vectors must be of the same size");
    }
    if constexpr (std::is_same_v<T, float>){
      cblas_saxpy(data_.size(), alpha, v.data_.data(), 1, data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_daxpy(data_.size(), alpha, v.data_.data(), 1, data_.data(), 1);
    }    			
  }


  Vector& operator += (const Vector<T>& v){
    // Vector addition: this += v
    
    if (data_.size() != v.data_.size()) {
      throw std::invalid_argument("Vectors must be of the same size");
    }
    add(1.0, v);
    return *this;
  }


  Vector& operator -= (const Vector<T>& v){
    // Vector addition: this -= v    
    if (data_.size() != v.data_.size()) {
      throw std::invalid_argument("Vectors must be of the same size");
    }
    add(-1.0, v);
    return *this;
  }
  
  void subtract(const Vector<T>& v) {
    // Vector subtraction: this = this - v    
    if (data_.size() != v.data_.size()) {
      throw std::invalid_argument("Vectors must be of the same size");
    }
    if constexpr (std::is_same_v<T, float>){
      cblas_saxpy(data_.size(), -1.0, v.data_.data(), 1, data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      cblas_daxpy(data_.size(), -1.0, v.data_.data(), 1, data_.data(), 1);
    }
  }

  T dot(const Vector<T>& v) const {
    // Dot product with another vector    
    if (data_.size() != v.data_.size()) {
      throw std::invalid_argument("Vectors must be of the same size");
    }
    if constexpr (std::is_same_v<T, float>){
      return cblas_sdot(data_.size(), data_.data(), 1, v.data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      return cblas_ddot(data_.size(), data_.data(), 1, v.data_.data(), 1);	
    }
  }

public:////////////////////////////////////////////////  getter and setter

  T* data(){
    return data_.data();
  }

public://///////////////////////////////////////////////// Memory management
  
  
public://////////////////////////////////////////// Iterators and Range-based Access

  //typedef T* iterator;
  using iterator = T*;
  //typedef const T* const_iterator;
  using const_iterator = const T*;
  

  iterator begin (){
    return data_.data();
  };

  /**
   * Return constant iterator to the start of
   * the vectors.
   */
  const_iterator begin () const{
    return data_.data();
  };

  /**
   * Return an iterator pointing to the
   * element past the end of the array.
   */
  iterator end (){
    return data_.data() + data_.size();
  };

  /**
   * Return a constant iterator pointing to
   * the element past the end of the array.
   */
  const_iterator end () const{
    return data_.data() + data_.size();
  };
  //@}

  

  
public://////////////////////////////////////////// Utility of Query functions
  /// get size of the vector
  size_t size() const {
    return data_.size();
  }

  // L2 Norm of the vector
  T l2_norm() const {
    if constexpr (std::is_same_v<T, float>){
      return cblas_snrm2(data_.size(), data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      return cblas_dnrm2(data_.size(), data_.data(), 1);
    }
  }

  // L1 Norm of the vector
  T l1_norm() const {
    if constexpr (std::is_same_v<T, float>){
      return cblas_sasum(data_.size(), data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      return cblas_dasum(data_.size(), data_.data(), 1);
    }
  }

  // L8 Norm of the vector
  T l_inf_norm() const {
    if constexpr (std::is_same_v<T, float>){
      return cblas_isamax(data_.size(), data_.data(), 1);
    }
    else if constexpr (std::is_same_v<T, double>){
      return cblas_idamax(data_.size(), data_.data(), 1);
    }
  }
  
  

  
public://///////////////////////////////////////////////////////// Interoperability and IO
  // Print vector elements
  void print() const {
    for (const auto& elem : data_) {
      std::cout << elem << " ";
    }
    std::cout << std::endl;
  }

  // Block write function
  void block_write(std::ostream &out) const {
    if (!out.good()) {
      throw std::runtime_error("Output stream is not good.");
    }

    size_t size = data_.size();
    
    // Write the size of the vector
    out.write(reinterpret_cast<const char*>(&size), sizeof(size));
    if (!out.good()) {
      throw std::runtime_error("Failed to write vector size to output stream.");
    }

    // Write the elements of the vector
    out.write(reinterpret_cast<const char*>(data_.data()), size * sizeof(T));
    if (!out.good()) {
      throw std::runtime_error("Failed to write vector data to output stream.");
    }
  }

  // Block read function
  void block_read(std::istream &in) {
    if (!in.good()) {
      throw std::runtime_error("Input stream is not good.");
    }

    size_t size;
    
    // Read the size of the vector
    in.read(reinterpret_cast<char*>(&size), sizeof(size));
    if (!in.good() || in.gcount() != sizeof(size)) {
      throw std::runtime_error("Failed to read vector size from input stream.");
    }

    // Resize the internal data vector
    data_.resize(size);

    // Read the elements of the vector
    in.read(reinterpret_cast<char*>(data_.data()), size * sizeof(T));
    if (!in.good() || in.gcount() != static_cast<std::streamsize>(size * sizeof(T))) {
      throw std::runtime_error("Failed to read vector data from input stream.");
    }
  }
  

  ////////////////////////////////////////////////////////////////////////////////////////////

private:
  std::vector<T> data_;
};
#endif//_VECTOR_H_
