
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 2.1 release candidate  (release date: 31 december 2016)
 
 Copyright (c), 2016 - GEOtop Foundation
 
 This file is part of GEOtop 2.1 
 
 GEOtop 2.1  is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 2.1  is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <vector>
#include <string>
#define __MATHOPTIM_H__
#include <meteoio/MeteoIO.h>

#if defined(__linux) && !defined(ANDROID) && !defined(__CYGWIN__)
#include <execinfo.h>
/* #define ENABLE_PRINT_STACK_FRAME 1 */

#ifdef ENABLE_PRINT_STACK_FRAME
#define PRINT_FRAME_STACK \
void *__array[2]; \
size_t __size; \
char **__strings; \
\
__size = backtrace (__array, 2); \
__strings = backtrace_symbols (__array, __size); \
\
printf ("Printing mangled stack frames: %zd\n", __size); \
\
for (size_t __i = 0; __i < __size; __i++) \
printf ("%s\n", __strings[__i]); \
\
free (__strings);
#else
#define PRINT_FRAME_STACK
#endif
#endif

#if defined _WIN32 || defined __MINGW32__
#include <windows.h>
#endif

template <class T> class GeoVector {
	public:
		GeoVector(const size_t& asize=0);

		/**
		* A constructor that creates a vector filled with constant values
		* @param asize size of the new array
		* @param init initial value to fill the vector with
		*/
		GeoVector(const size_t& asize, const T& init);

		T& operator [](const size_t& index);
		const T operator [](const size_t& index) const;
		T& at(const size_t& index);
		const T at(const size_t& index) const;

		size_t size() const;
		void resize(const size_t& asize);
		void resize(const size_t& asize, const T& init);
		void clear();

		/**
		* Calling this void procedure sets all elements in GeoVector to val
		* (the size is not affected)
		* @param val Value of type T 
		*/
		void reset(const T& val);

		/**
		* Calling this void procedure sets all elements in GeoVector to val
		* except the ones that already have value val_to_omit (the size is not affected)
		* @param val Value of type T that will be copied to all elements, except if
		* @param val_to_omit Value of type T that will be preserved
		*/
		void reset(const T& val, const T& val_to_omit);

		void setMetaData(const std::string& name="unknown", const std::string& unit="unknown", const size_t& precision=3);

		std::vector<T> data;
		std::string name;
		std::string unit;
		size_t precision;
};

template<class T> GeoVector<T>::GeoVector(const size_t& asize)
{
	resize(asize);
	setMetaData();
}

template<class T> GeoVector<T>::GeoVector(const size_t& asize, const T& init)
{
	resize(asize, init);
	setMetaData();
}

template<class T> void GeoVector<T>::setMetaData(const std::string& iname, const std::string& iunit, const size_t& iprecision)
{
	name = iname;
	unit = iunit;
	precision = iprecision;
}

template<class T> void GeoVector<T>::reset(const T& val)
{
  const size_t data_size = data.size();
	for (size_t ii=0; ii<data_size; ii++) {
		data[ii] = val;
	}
}

template<class T> void GeoVector<T>::reset(const T& val, const T& val_to_omit)
{
	const double data_size = data.size();
	for (size_t ii=0; ii<data_size; ii++) {
		if (data[ii] != val_to_omit) {
			data[ii] = val;
		}
	}
}


template<class T> size_t GeoVector<T>::size() const
{
#if defined(__linux)
    PRINT_FRAME_STACK
#endif
	return data.size();
}

template<class T> void GeoVector<T>::resize(const size_t& asize)
{
	data.resize(asize);
}

template<class T> void GeoVector<T>::resize(const size_t& asize, const T& init)
{
	resize(asize);
	std::fill(data.begin(), data.end(), init);
}

template<class T> inline T& GeoVector<T>::at(const size_t& index)
{
	return data.at(index);
}

template<class T> inline const T GeoVector<T>::at(const size_t& index) const
{
	return data.at(index);
}

template<class T> inline T& GeoVector<T>::operator [](const size_t& index)
{
	return data[index];
}

template<class T> inline const T GeoVector<T>::operator [](const size_t& index) const
{
	return data[index];
}

template<class T> void GeoVector<T>::clear()
{
	resize(0);
}


template <class T> class GeoMatrix : public mio::Array2D<T> {

	public:

		GeoMatrix();

		/**
		* A constructor that creates a matrix of a given size
		* @param anx number of columns of the new matrix
		* @param any number of rows of the new matrix
		*/
    GeoMatrix(const size_t& anx, const size_t &any);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param init initial value to fill the array with
		*/
    GeoMatrix(const size_t& anx, const size_t& any, const T& init);

    size_t getRows() const;
    size_t getCols() const;

		/**
		* Calling this void procedure sets all elements in GeoMatrix to val
		* the size is not affected
		* @param val Value of type T 
		*/
		void reset(const T& val);

		/**
		* Calling this void procedure sets all elements in GeoMatrix to val
		* except the ones that already have value val_to_omit (the size is not affected)
		* @param val Value of type T that will be copied to all elements, except if
		* @param val_to_omit Value of type T that will be preserved
		*/
		void reset(const T& val, const T& val_to_omit);

		void setMetaData(const std::string& name="unknown", const std::string& unit="unknown", const size_t& precision=3);

		std::string name;
		std::string unit;
		size_t precision;
};

template<class T> GeoMatrix<T>::GeoMatrix() : mio::Array2D<T>()
{
	setMetaData();
}

template<class T> GeoMatrix<T>::GeoMatrix(const size_t& anx, const size_t &any) : mio::Array2D<T>(anx, any)
{
	setMetaData();
}

template<class T> GeoMatrix<T>::GeoMatrix(const size_t &anx, const size_t& any, const T& init)  : mio::Array2D<T>(anx, any, init)
{
	setMetaData();
}

template<class T> void GeoMatrix<T>::setMetaData(const std::string& iname, const std::string& iunit, const size_t& iprecision)
{
	name = iname;
	unit = iunit;
	precision = iprecision;
}

template<class T> void GeoMatrix<T>::reset(const T& val)
{
	for (unsigned int ii=0; ii<this->nx; ii++) {
		for (unsigned int jj=0; jj<this->ny; jj++) {
			this->operator()(ii,jj) = val;
		}
	}
}

template<class T> void GeoMatrix<T>::reset(const T& val, const T& val_to_omit)
{
	for (unsigned int ii=0; ii<this->nx; ii++) {
		for (unsigned int jj=0; jj<this->ny; jj++) {
			if (this->operator()(ii,jj) != val_to_omit) {
				this->operator()(ii,jj) = val;
			}
		}
	}
}

template<class T> size_t GeoMatrix<T>::getRows() const
{
#if defined(__linux)
    PRINT_FRAME_STACK
#endif
	return this->nx;
}

template<class T> size_t GeoMatrix<T>::getCols() const
{
#if defined(__linux)
    PRINT_FRAME_STACK
#endif
	return this->ny;
}


template <class T> class GeoTensor : public mio::Array3D<T> {

	public:

		GeoTensor();

		/**
		* A constructor that creates a 3d array of a given size
		* @param anx number of columns of the new 3d array
		* @param any number of rows of the new 3d array
		* @param anz number of depths of the new 3d array
		*/
    GeoTensor(const size_t& anx, const size_t& any, const size_t& anz);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param anz number of depths of the new 3d array
		* @param init initial value to fill the array with
		*/
    GeoTensor(const size_t& anx, const size_t& any, const size_t& anz, const T& init);

    size_t getDh() const;
    size_t getCh() const;
    size_t getRh() const;

		/**
		* Calling this void procedure sets all elements in GeoTensor to val
		* the size is not affected
		* @param val Value of type T 
		*/
		void reset(const T& val);

		/**
		* Calling this void procedure sets all elements in GeoTensor to val
		* except the ones that already have value val_to_omit (the size is not affected)
		* @param val Value of type T that will be copied to all elements, except if
		* @param val_to_omit Value of type T that will be preserved
		*/
		void reset(const T& val, const T& val_to_omit);


		void setMetaData(const std::string& name="unknown", const std::string& unit="unknown", const size_t& precision=3);

		std::string name;
		std::string unit;
		size_t precision;
};

template<class T> GeoTensor<T>::GeoTensor() : mio::Array3D<T>()
{
	setMetaData();
}

template<class T> GeoTensor<T>::GeoTensor(const size_t &anx, const size_t &any, const size_t &anz) : mio::Array3D<T>(anx, any, anz)
{
	setMetaData();
}

template<class T> GeoTensor<T>::GeoTensor(const size_t &anx, const size_t &any, const size_t &anz, const T& init)  : mio::Array3D<T>(anx, any, anz, init)
{
	setMetaData();
}

template<class T> void GeoTensor<T>::setMetaData(const std::string& iname, const std::string& iunit, const size_t& iprecision)
{
	name = iname;
	unit = iunit;
	precision = iprecision;
}

template<class T> void GeoTensor<T>::reset(const T& val)
{
	for (unsigned int ii=0; ii<this->nx; ii++) {
		for (unsigned int jj=0; jj<this->ny; jj++) {
			for (unsigned int kk=0; kk<this->nz; kk++) {
				this->operator()(ii,jj,kk) = val;
			}
		}
	}
}

template<class T> void GeoTensor<T>::reset(const T& val, const T& val_to_omit)
{
	for (unsigned int ii=0; ii<this->nx; ii++) {
		for (unsigned int jj=0; jj<this->ny; jj++) {
			for (unsigned int kk=0; kk<this->nz; kk++) {
				if (this->operator()(ii,jj,kk) != val_to_omit) {
					this->operator()(ii,jj,kk) = val;
				}
			}

		}
	}
}

template<class T> size_t GeoTensor<T>::getDh() const
{
#if defined(__linux)
    PRINT_FRAME_STACK
#endif
	return this->nx;
}

template<class T> size_t GeoTensor<T>::getRh() const
{
#if defined(__linux)
    PRINT_FRAME_STACK
#endif
	return this->ny;
}

template<class T> size_t GeoTensor<T>::getCh() const
{
#if defined(__linux)
    PRINT_FRAME_STACK
#endif
	return this->nz;
}

class TInit {  	   /*header of maps in fluid turtle format*/
public:
     GeoVector<double> U;  /*dx,dy*/
    GeoVector<double> V;  /*sign of novalue,novalue*/
};



#endif
