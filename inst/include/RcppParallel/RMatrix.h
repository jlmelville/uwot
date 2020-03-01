#ifndef __RCPP_PARALLEL_RMATRIX__
#define __RCPP_PARALLEL_RMATRIX__

#include <cstddef>
#include <iterator>

namespace RcppParallel {

template <typename T>
class RMatrix {
public:
   class Row {
   
   public:   
      
      template <typename V>
      class row_iterator 
         : public std::iterator<std::random_access_iterator_tag, V, std::size_t> {
      
      public:
         inline row_iterator(Row& row, std::size_t i)
            : start_(row.start_), parentNrow_(row.parent_.nrow()), index_(i)
         {
         }
         
         inline row_iterator(V* start, std::size_t parentNrow, std::size_t index)
            : start_(start), parentNrow_(parentNrow), index_(index)
         {      
         }
         
         inline row_iterator(const row_iterator& other) 
            : start_(other.start_), 
              parentNrow_(other.parentNrow_), 
              index_(other.index_)
         {
         }
         
         inline row_iterator& operator++() { 
            index_++;
            return *this;
         }
         
         inline row_iterator operator++(int) {
            row_iterator tmp(*this); 
            operator++(); 
            return tmp;
         }
         
         inline row_iterator& operator--() {
            index_-- ;
            return *this ;
         }
         
         inline row_iterator operator--(int) {
            row_iterator tmp(*this);
            index_-- ;
            return tmp ;
         }

         row_iterator operator+(std::size_t n) const { 
            return row_iterator(start_, parentNrow_ ,index_ + n ) ; 
         }
         row_iterator operator-(std::size_t n) const { 
            return row_iterator(start_, parentNrow_, index_ - n ) ; 
         }
         
         std::size_t operator+(const row_iterator& other) const {
            return index_ + other.index_;
         }
         
         std::size_t operator-(const row_iterator& other) const { 
            return index_ - other.index_ ; 
         }
         
         row_iterator& operator+=(std::size_t n) { index_ += n ; return *this; }
         row_iterator& operator-=(std::size_t n) { index_ -= n ; return *this; }
         
         bool operator==(const row_iterator& other) const { return index_ == other.index_; }
         bool operator!=(const row_iterator& other) const { return index_ != other.index_; }
         bool operator<(const row_iterator& other) const { return index_ < other.index_; }
         bool operator>(const row_iterator& other) const { return index_ > other.index_; }
         bool operator<=(const row_iterator& other) const { return index_ <= other.index_; }
         bool operator>=(const row_iterator& other) const { return index_ >= other.index_; }
         

         inline V& operator*() { return start_[index_ * parentNrow_]; }
         
         inline V* operator->() { return &(start_[index_ * parentNrow_]); }
      
         inline V& operator[](int i) { return start_[(index_+i) * parentNrow_]; }
         
      private:
         V* start_;
         std::size_t parentNrow_;
         std::size_t index_;
      };
      
      typedef row_iterator<T> iterator;
      typedef row_iterator<const T> const_iterator;
      
      inline Row(RMatrix& parent, std::size_t i)
         : parent_(parent),
           start_(parent.begin() + i)
      {
      }
      
      inline Row(const Row& other)
         : parent_(other.parent_),
           start_(other.start_)
      {        
      }
      
      inline iterator begin() {
         return iterator(*this, 0);
      }
      
      inline iterator end() {
         return iterator(*this, parent_.ncol());
      }
      
      inline const_iterator begin() const {
         return const_iterator(*this, 0);
      }
      
      inline const_iterator end() const {
         return const_iterator(*this, parent_.ncol());
      }
      
      inline size_t length() const {
         return parent_.ncol();
      }

      inline size_t size() const {
        return parent_.ncol();
      }

      inline T& operator[](std::size_t i) {
        return start_[i * parent_.nrow()];
      }
      
      inline const T& operator[](std::size_t i) const {
        return start_[i * parent_.nrow()];
      }
              
   private:
      RMatrix& parent_;
      T* start_;
   };
   
   class Column {
   
   public:
   
      typedef T* iterator;
      typedef const T* const_iterator;
   
      inline Column(RMatrix& parent, std::size_t i) 
         : begin_(parent.begin() + (i * parent.nrow())),
           end_(begin_ + parent.nrow())
      {   
      }
      
      inline Column(const Column& other) 
         : begin_(other.begin_), end_(other.end_)
      {   
      }
      
      inline Column& operator=(const Column& rhs) {
         begin_ = rhs.begin_;
         end_ = rhs.end_;
         return *this;
      }
      
      inline iterator begin() { return begin_; }
      inline iterator end() { return end_; }
      
      inline const_iterator begin() const { return begin_; }
      inline const_iterator end() const { return end_; }
      
      inline size_t length() const { return end_ - begin_; }
      inline size_t size() const { return end_ - begin_; }
      
      inline T& operator[](std::size_t i) {
        return *(begin_ + i);
      }
      
      inline const T& operator[](std::size_t i) const {
        return *(begin_ + i);
      }
      
   private:
      T* begin_;
      T* end_;
   };

   typedef T* iterator;
   typedef const T* const_iterator;

   template <typename Source>
   inline explicit RMatrix(const Source& source) 
      : data_(const_cast<Source&>(source).begin()),
        nrow_(source.nrow()),
        ncol_(source.ncol())
   {
   }

   inline RMatrix(const T* data, std::size_t nrow, std::size_t ncol) 
      : data_(data), nrow_(nrow), ncol_(ncol) 
   {
   }
   
   inline iterator begin() { return data_; }
   inline iterator end() { return data_ + length(); }
   
   inline const_iterator begin() const { return data_; }
   inline const_iterator end() const { return data_ + length(); }
     
   inline std::size_t length() const { return nrow_ * ncol_; }  
     
   inline std::size_t nrow() const { return nrow_; }
   inline std::size_t ncol() const { return ncol_; }
   
   inline T& operator()(std::size_t i, std::size_t j) {
      return *(data_ + (i + j * nrow_));
   }
   
   inline const T& operator()(std::size_t i, std::size_t j) const {
      return *(data_ + (i + j * nrow_));
   }
   
   inline Row row(std::size_t i) {
      return Row(*this, i);
   }
   
   inline const Row row(std::size_t i) const {
      return Row(*const_cast<RMatrix*>(this), i);
   }
   
   inline Column column(std::size_t i) {
      return Column(*this, i);
   }
   
   inline const Column column(std::size_t i) const {
      return Column(*const_cast<RMatrix*>(this), i);
   }
   
   inline T& operator[](std::size_t i) {
      return *(data_ + i);
   }
   
   inline const T& operator[](std::size_t i) const {
      return *(data_ + i);
   }
   
private:
   T* data_;
   std::size_t nrow_;
   std::size_t ncol_;
};

} // namespace RcppParallel

#endif // __RCPP_PARALLEL_RMATRIX__
