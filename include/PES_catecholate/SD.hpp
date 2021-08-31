# if ! defined(_SD_H)
# define _SD_H

# include <PES_catecholate/constants.hpp>
# include <iostream>
# include <cmath>
# include <vector>

/** function template to calculate norm of a vector **/

template<typename T>
inline double Vec_Norm(const std::vector<T>& vect)	{
	double norm = pow(vect[0],2);
	for (unsigned int i=1; i<vect.size(); i++)	{
		norm += pow(vect[i],2);
	}
	norm = sqrt(norm);
	return norm;
}

/** function template to make all elements zero **/

template<typename T>
void make_vect_zero(std::vector<T>& vect)	{
	for (size_t i=0; i<vect.size(); i++)	{
		vect[i] = 0.;
	}
}

/** function template to multiply all elements of vector with a scalar **/

template<typename T>
inline void multiply_w_const(std::vector<T>& inputvect,T val)	{
	for (size_t i=0; i<inputvect.size(); i++)	{
		inputvect[i] *= val;
	}
}

/** function template to add second vector to first vector in the argument **/

template<typename T>
inline void add_2_vect(std::vector<T>& vect1,const std::vector<T>& vect2)	{
	for (size_t i=0; i<vect1.size(); i++)	{
		vect1[i] += vect2[i];
	}
}

/*** function template to print out 1D vector on terminal ***/

template<typename T>
inline void printonedimvect(const std::vector<T>& inputvect)	{
	for (size_t i=0; i<inputvect.size(); i++)	{
		std::cout	<< std::scientific	<<	inputvect[i]	<< "\n";
	}
}

/*** function template to print out 2D vector on terminal ***/

template<typename T>
inline void printtwodimvect(const std::vector<std::vector<T>>& inputvect)	{
	for (size_t i=0; i<inputvect.size(); i++)	{
		for (size_t j=0; j<inputvect[i].size(); j++)	{
			std::cout	<< std::scientific	<< inputvect[i][j]	<< "	";
		}
		std::cout	<< "\n";
	}
}

/** function to split a std::vector into 2D std::vector<std::vector<T>> 
 * with 1D std::vector and size of each partition as inputs **/

template<typename T>
inline std::vector<std::vector<T>> split_onedim_vector(const std::vector<T>& onedimvect,int size)	{
	std::vector<std::vector<T>> rtrn_twodim_vect;
	std::vector<T> tmpvect;
	int counter = 0;
	for (size_t i=0; i<onedimvect.size(); i++)	{
		counter++;
		tmpvect.emplace_back(onedimvect[i]);
		if (counter == size)	{
			rtrn_twodim_vect.push_back(tmpvect);
			counter = 0;
			tmpvect.clear();
		}
	}
	return rtrn_twodim_vect;
}

/** function template to calculate first derivative using 5-pt numerical differentiation **/

template<typename T>
std::vector<double> get_first_drv(const std::vector<std::vector<T>>& funcvect,T stepsize)	{
	std::vector<double> first_drv_vect(funcvect.size());
	for (size_t i=0; i<first_drv_vect.size(); i++)	{
		if (funcvect[i].size() != 5)	{
			throw "No of data points not suitable for 5-Pt differentiation!";
			exit(1);
		}
		first_drv_vect[i] = static_cast<double>((1/stepsize))*(SD_PARAMS::first_drv_coeff_1*funcvect[i][0]+
				SD_PARAMS::first_drv_coeff_2*funcvect[i][1]+SD_PARAMS::first_drv_coeff_3*funcvect[i][2]+
				SD_PARAMS::first_drv_coeff_4*funcvect[i][3]+SD_PARAMS::first_drv_coeff_5*funcvect[i][4]); 
	}
	return first_drv_vect;
}


/** function for generating displacement 
 * grids for all coordinates for gradient calculation 
 *  using 5-Pt numerical differentiation **/

std::vector<std::vector<double>> get_all_onedim_grids_4_grad(const std::vector<double>& curgrid);

/** function to generate multidimensional grids 
 * for function evaluation at displaced grids 
 * for all coordinates to be used for gradient 
 * calculation with 5-Pt numerical differentiation ****/

std::vector<std::vector<double>> get_multidim_grid_from_all_onedim_grids(const std::vector<std::vector<double>>& all_onedim_grids,const std::vector<double>& const_coord_vals);



# endif
