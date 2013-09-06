
/* STATEMENT:
 
 GEOtop MODELS THE ENERGY AND WATER FLUXES AT THE LAND SURFACE
 GEOtop 1.225 'Moab' - 9 Mar 2012
 
 Copyright (c), 2012 - Stefano Endrizzi
 
 This file is part of GEOtop 1.225 'Moab'
 
 GEOtop 1.225 'Moab' is a free software and is distributed under GNU General Public License v. 3.0 <http://www.gnu.org/licenses/>
 WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
 
 GEOtop 1.225 'Moab' is distributed as a free software in the hope to create and support a community of developers and users that constructively interact.
 If you just use the code, please give feedback to the authors and the community.
 Any way you use the model, may be the most trivial one, is significantly helpful for the future development of the GEOtop model. Any feedback will be highly appreciated.
 
 If you have satisfactorily used the code, please acknowledge the authors.
 
 */
#ifndef DATASTRUCTS_H
#define DATASTRUCTS_H

#include <vector>
#include <string>
#include <meteoio/MeteoIO.h>

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
	const double data_size = data.size();
	for (unsigned int ii=0; ii<data_size; ii++) {
		data[ii] = val;
	}
}

template<class T> void GeoVector<T>::reset(const T& val, const T& val_to_omit)
{
	const double data_size = data.size();
	for (unsigned int ii=0; ii<data_size; ii++) {
		if (data[ii] != val_to_omit) {
			data[ii] = val;
		}
	}
}


template<class T> size_t GeoVector<T>::size() const
{
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
		GeoMatrix(const unsigned int& anx, const unsigned int& any);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param init initial value to fill the array with
		*/
		GeoMatrix(const unsigned int& anx, const unsigned int& any, const T& init);

		unsigned int getRows() const;
		unsigned int getCols() const;

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

template<class T> GeoMatrix<T>::GeoMatrix(const unsigned int& anx, const unsigned int& any) : mio::Array2D<T>(anx, any)
{
	setMetaData();
}

template<class T> GeoMatrix<T>::GeoMatrix(const unsigned int& anx, const unsigned int& any, const T& init)  : mio::Array2D<T>(anx, any, init)
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

template<class T> unsigned int GeoMatrix<T>::getRows() const 
{
	return this->nx;
}

template<class T> unsigned int GeoMatrix<T>::getCols() const 
{
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
		GeoTensor(const unsigned int& anx, const unsigned int& any, const unsigned int& anz);

		/**
		* A constructor that creates an array filled with constant values
		* @param anx number of columns of the new array
		* @param any number of rows of the new array
		* @param anz number of depths of the new 3d array
		* @param init initial value to fill the array with
		*/
		GeoTensor(const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init);

		unsigned int getDh() const;
		unsigned int getCh() const;
		unsigned int getRh() const;

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

template<class T> GeoTensor<T>::GeoTensor(const unsigned int& anx, const unsigned int& any, const unsigned int& anz) : mio::Array3D<T>(anx, any, anz)
{
	setMetaData();
}

template<class T> GeoTensor<T>::GeoTensor(const unsigned int& anx, const unsigned int& any, const unsigned int& anz, const T& init)  : mio::Array3D<T>(anx, any, anz, init)
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

template<class T> unsigned int GeoTensor<T>::getDh() const 
{
	return this->nx;
}

template<class T> unsigned int GeoTensor<T>::getRh() const 
{
	return this->ny;
}

template<class T> unsigned int GeoTensor<T>::getCh() const 
{
	return this->nz;
}

#endif
