#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <algorithm>
#include <iostream>
#include <cstring>
#include <string>
#include <exception>
using std::fill;
using std::memcpy;
using std::calloc;
using std::free;
using std::cin;
using std::cout;
using std::to_string;

template <class T>
class matrix{
    public:
    int nrow;
    int ncol;
    T* data;

    // Constructors
    matrix ();
    matrix (size_t nRow, size_t nCol);
    matrix (size_t nRow, size_t nCol, const T& e);
    ~matrix ();

    // Copy constructors
    matrix (const matrix<T>& m);

    // Assignment operators
    matrix<T>& operator= (const T& e);
    matrix<T>& operator= (const matrix<T>& mat);

    size_t size (){ return nrow*ncol; }

    T& operator() (int i, int j) const;
    T& operator[] (int i) const;


};




template <class T> inline
matrix<T>::matrix (){
    nrow=0;
    ncol=0;
    data = nullptr;
}

template <class T> inline
matrix<T>::matrix (size_t nRow, size_t nCol){
    nrow=nRow;
    ncol=nCol;
    data =  (T*)calloc(nrow*ncol,sizeof(T) );
    fill(data,data+nrow*ncol,0);
}

template <class T> inline
matrix<T>::matrix (size_t nRow, size_t nCol, const T& e){
    nrow=nRow;
    ncol=nCol;
    data =  (T*)calloc(nrow*ncol,sizeof(T) );
    fill(data,data+nrow*ncol,e);
}


template <class T> inline
matrix<T>::matrix (const matrix<T>& m){
    nrow=m.nrow;
    ncol=m.ncol;
    data =  (T*)calloc(nrow*ncol,sizeof(T) );
    memcpy(data,m.data,nrow*ncol*sizeof(T) );
}


template <class T> inline
matrix<T>::~matrix (){
    nrow=0.0;
    ncol=0.0;
    if(data){
        free(data);
        data=nullptr;
    }
}


template <class T>
matrix<T>& matrix<T>::operator= (const T& el){
    fill(data,data+nrow*ncol,el);
    return *this;
}


template <class T>
matrix<T>& matrix<T>::operator= (const matrix<T>& mat){
    nrow=mat.nrow;
    ncol=mat.ncol;
    data=(T*)calloc(nrow*ncol,sizeof(T));
    memcpy(data,mat.data,nrow*ncol*sizeof(double));
    // fill(data,data+nrow*ncol,el);
    return *this;
}


template <class T> inline T&
matrix<T>::operator() (int i, int j) const{
    int row=i, column=j;
    if(i<-nrow || i>nrow-1 || j<-ncol || j>ncol-1) throw std::out_of_range("Index out of range in matrix::matrix("+to_string(i)+","+to_string(j)+")!");
    if(i<0) row=nrow+i;
    if(j<0) column=ncol+j;
    return data[row*ncol+column ];
}


template <class T> inline T&
matrix<T>::operator[] (int i) const{
    int row=i;
    if(i<-nrow || i>nrow-1) throw std::out_of_range("Index out of range in matrix::matrix["+to_string(i)+"]!");
    if(i<0) row=nrow+i;
    return data[row*ncol ];
}

#endif

