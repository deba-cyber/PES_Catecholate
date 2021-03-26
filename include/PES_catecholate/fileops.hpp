# if ! defined (FILEOPS)
# define FILEOPS

# include <iostream>
# include <fstream>
# include <sstream>
# include <string>
# include <vector>
# include <algorithm>



/****************************************************
 * Function template to free up memory allocated 
 * for vector of any data type	********************
 * **************************************************/

template<typename T>
inline auto empty_swap(std::vector<T> & vec)	{
	std::vector<T>().swap(vec);
	return vec;
}

/****************************************************
 * Function template to return specific queried row *
 * or column elements from passed one dimensional 
 * vector (which is obtained from two dimensional 
 * file data ........................***************
 * **************************************************/

template<typename T>
inline auto get_specific_row_or_col_elements(std::vector<T> & vec_onedim,
		int N_row_or_col, int N_tot_rows, int N_tot_cols, 
		const std::string s1)	
{
	/** Indexing for N_row_or_col is 0-based i.e. if 
	 *	N_row_or_col is 100 , actually 101-th row is desired **/
	std::vector<T> rtrn_sp_row_or_col_elems;
	std::vector<int> col_iterator_vect;		// index vector over which loop will run for column
	int start_cnt;
	int end_cnt;
	int cur_row = 0;
	int cur_id_4_col = N_row_or_col;		// index for iterating over for extracting column 
	if (s1 == "ROW")	{
		start_cnt = N_row_or_col*N_tot_cols;
		end_cnt = (N_row_or_col+1)*N_tot_cols;
		for (int i=start_cnt; i<end_cnt; i++)	{
			rtrn_sp_row_or_col_elems.emplace_back(vec_onedim[i]);
		}
		return rtrn_sp_row_or_col_elems;
	}
	else if (s1 == "COLUMN")	{
		while (cur_row != N_tot_rows)	{
			rtrn_sp_row_or_col_elems.emplace_back(vec_onedim[cur_id_4_col]);
			cur_id_4_col += N_tot_cols;
			cur_row +=1;
		}
		return rtrn_sp_row_or_col_elems;
	}
}

/*******************************************************
 * Function template for extracting all lines from
 * file and storing in a multidimensional vector
 * Useful when each line in the file has different 
 * number of elements ..
 * lineno_vect has line numbers to be extracted (1-based
 * indexing)
 * each element of the multidimensional vector has length
 * equal to number of elements in the line in the file
 * passing dummy vector to keep template structure
*******************************************************/

template<typename T>
inline std::vector<std::vector<T>> Get_all_lines_multvect(const std::vector<int>& lineno_vect,
		const std::vector<T>& dummyvec,const std::string& fstring)	{
	std::string fileext = ".dat";
	std::string line;
	std::string filepath = fstring + fileext; 
	std::vector<std::vector<T>> rtrn_multline_multvect;
	std::ifstream myfile(filepath);
	int l_cnt = 0;
	for (size_t i=0; i<lineno_vect.size(); i++)	{
		int cur_lineno = lineno_vect[i];
		std::vector<T> tmp_vect;
		if (myfile.is_open())	{
			while (l_cnt != cur_lineno && std::getline(myfile,line))	{
				++ l_cnt;
			}
			if (l_cnt == cur_lineno)	{
				std::stringstream lstream(line);
				T val;
				while (lstream >> val)	{
					tmp_vect.emplace_back(val);
				}
				rtrn_multline_multvect.push_back(tmp_vect);
			}
		}
	}
	return rtrn_multline_multvect;
}


/*******************************************************
 * Function template to subtract a vector from another
********************************************************/

template<typename T>
inline std::vector<T> subtract_vect(const std::vector<T>& vect1,
		const std::vector<T>& vect2)	{
	std::vector<T> rtrn_vect(vect1.size());
	if (vect1.size() != vect2.size())	{
		throw "Vector sizes are not equal!";
	}
	else {
		std::vector<T> rtrn_vect(vect1.size());
		for (size_t i=0; i<rtrn_vect.size(); i++)	{
			rtrn_vect[i] = vect1[i]-vect2[i];
		}
	return rtrn_vect;
	}
}


/*******************************************************
Function template to calculate norm square of a vector
*******************************************************/

template<typename T>
inline T get_norm_sqr(const std::vector<T>& vect)	{
	T norm_sqr = 0.;
	for (size_t i=0; i<vect.size(); i++)	{
		norm_sqr += vect[i]*vect[i];
	}
	return norm_sqr;
}

/*********************************************************
Function template to calculate dot product of two vectors
*********************************************************/

template<typename T>
inline T get_dot_pdt(const std::vector<T>& vect1,
		const std::vector<T>& vect2)	{
	T dot_pdt = 0.;
	if (vect1.size() != vect2.size())	{
		throw "Vector sizes are not equal!";
		exit(1);
	}
	else {
		for (size_t i=0; i<vect1.size(); i++)	{
			dot_pdt += vect1[i]*vect2[i];
		}
	}
	return dot_pdt;
}

/*******************************************************
 * Function template for generating transpose of a matrix
 * getting input as vector and returning as vector
*******************************************************/

template<typename T>
inline std::vector<T> get_mat_transpose_as_vect(const std::vector<T>& input_mat,
		int N_row,int N_col)	{
	std::vector<T> rtrn_transpose;
	rtrn_transpose.reserve(input_mat.size());
	for (int i=0; i<N_col; i++)	{
		for (int j=0; j<N_row; j++)	{
			rtrn_transpose.emplace_back(input_mat[i+j*N_col]);
		}
	}
	return rtrn_transpose;
}

/*******************************************************
 * Function template to return matrix-vector product ***
 * where the matrix is passed in argument as 1-dim 
 * vector	by reference .............................
 * No of rows and columns are also passed by values 
 * *****************************************************/

template<typename T>
inline std::vector<T> get_mat_vec_pdt(std::vector<T> & mat_as_vect, std::vector<T> & onedim_vect,
		int N_row, int N_col)
{
	std::vector<T> rtrn_matvec_pdt_vect;
	std::vector<T> tmp_process_vect;
	int onedim_vect_size = onedim_vect.size();
	T sum;
	if (N_col != onedim_vect_size)	{
		throw "Matrix-Vector multiplication not possible!";
	}
	else
	{
		for (int i=0; i<N_row; i++)	{
			tmp_process_vect = get_specific_row_or_col_elements(mat_as_vect,i,N_row,N_col,"ROW");
			sum = 0.;
			for (size_t j=0; j<tmp_process_vect.size(); j++)	{
				sum += tmp_process_vect[j]*onedim_vect[j];
			}
			rtrn_matvec_pdt_vect.emplace_back(sum);
		}
	}
	return rtrn_matvec_pdt_vect;
}


/************************************************
 * Template Function for extracting a specific 
 * line from a file of arbitrary length	.. ....
 * passing dummy vector having zero size of same 
 * data type for template structure ............ 
 * **********************************************/

template<typename T>
inline std::vector<T> get_specific_line(int lineno,const std::vector<T>& dummyvec,const std::string& s1)	{
	// lineno is given in 1-based indexing //
	// string s1 provides the filename without the extension //
	std::string line1;
	std::string ext = ".dat";
	std::vector<T> rtrn_sp_row;
	int l_cnt = 0;
	std::string path = s1 + ext;
	std::ifstream myfile(path);
	if (myfile.is_open())	{
		while (l_cnt< lineno && std::getline(myfile,line1))	{
			++l_cnt;
		}
		if (l_cnt == lineno)	{
			std::stringstream lstream(line1);
			T val;
			while (lstream >> val)	{
				rtrn_sp_row.emplace_back(val);
			}
		}
	}
	return rtrn_sp_row;
}

/**************************************************
 * Template Function for extracting specific lines
 * from a file	.................................
 * string s1 is passed by ref filename wo ext.
 * vector lineno_vect contains numbers in increasing
 * order corrsponsing to lines to be extracted	
 * lineno_vect has 1-based indexing	...............
 * returning as one dimensional vector ..........
 * passing dummy vector having zero size of same 
 * data type for template structure .............
 * ************************************************/

template<typename T>
inline std::vector<T> get_all_sp_lines(const std::vector<int> & lineno_vect, const std::vector<T>& dummyvec,
		const std::string& s1)		{
	std::string line1;
	std::string ext = ".dat";
	std::string path = s1 + ext;
	std::vector<T> rtrn_multiple_lines_vect;
	std::ifstream myfile(path);
	int cur_lineno;
	int l_cnt = 0;
	for (size_t i=0; i<lineno_vect.size(); i++)	{
		cur_lineno = lineno_vect[i];
		if (myfile.is_open())	{
			while (l_cnt!= cur_lineno	&& std::getline(myfile,line1))	{
				++l_cnt;
			}
			if (l_cnt == cur_lineno)	{
				std::stringstream lstream(line1);
				T val; 
				while (lstream >> val)	{
					rtrn_multiple_lines_vect.emplace_back(val);
				}
			}
		}
	}
	return rtrn_multiple_lines_vect;
}


/********************************************************************
 * function template for appending elements to an existing binary 
 * file ... 
 * data is passed by reference in an one dimensional vector
 * alongwith filename & mode which has to be "append" for this 
 * to work
 * ******************************************************************/

template<typename T>
inline void save_binary(std::vector<T> & vect,
		const std::string  filename, const std::string& mode)	{
	std::string ext = ".bin";
	if (mode == "append")	{
		std::ofstream outfile(filename+ext, std::ios::binary|std::ios::app);
		if (outfile.is_open())	{
			outfile.write(reinterpret_cast<char*>(&vect[0]),vect.size()*sizeof(T));
		}
		else	{
			throw "Error opening file for writing!\n";
			exit(1);
		}	
	}
}

/***********************************************************************
 * function template for extracting specific elements from a binary 
 * file ....
 * start_elem_no gives the element from which extraction starts
 * (1-based indexing)
 * vectsize gives number of elements to be extracted including the 
 * starting element ....................................................
 * *********************************************************************/


template<typename T>
inline std::vector<T> rtrn_vec_from_bin(std::vector<T> & dummyvec,
		const std::string& filename,int start_elem_no,
		int vectsize)	{
	/*****************************************
	 * passing a dummy vector of same data type 
	 * of zero size for keeping template structure
	 * ***************************************/
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext, std::ios::binary);
	std::vector<T> vect(vectsize);
	infile.seekg((start_elem_no-1)*sizeof(T),std::ios::beg);
	if (infile.is_open())	{
		infile.read(reinterpret_cast<char*>(&vect[0]),vectsize*sizeof(T));
	}
	else	{
		throw "Error opening file for writing!\n";
		exit(1);
	}
	return vect;
}

/***************************************************
 * function template to get number of elements in 
 * a binary file ..................................
 * *************************************************/


template<typename T>
inline int get_num_elems(std::vector<T> & dummyvec,
		const std::string& filename)	{
	/*****************************************
	 * passing a dummy vector of same data type 
	 * of size zero for keeping template structure
	 * ***************************************/
	std::ifstream infile;
	std::string ext = ".bin";
	infile.open(filename+ext,std::ios::binary);
	infile.seekg(0,std::ios::end);
	int filesize = infile.tellg();
	int num_elems = filesize/sizeof(T);
	return num_elems;
}


# endif
