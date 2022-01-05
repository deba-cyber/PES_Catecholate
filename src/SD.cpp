# include <PES_catecholate/SD.hpp>
# include <PES_catecholate/constants.hpp>
# include <iosfwd>
# include <vector>


STRIDE_ARR_BACK_N_FORTH::STRIDE_ARR_BACK_N_FORTH(const std::vector<int>& basis_size_vect)	{
	this->basis_size_vect = basis_size_vect;
}

std::vector<int> STRIDE_ARR_BACK_N_FORTH::stride_arr()	{
	/** preparing stride array **/
	int cur_product;
	int cur_index_pdt = 1;
	for (size_t j=1; j<basis_size_vect.size(); j++)	{
		cur_index_pdt*= basis_size_vect[j];
	}
	/** Initialising stride array with first element **/
	std::vector<int>stride_arr_init;
	stride_arr_init.emplace_back(cur_index_pdt);
	/** other elements of stride array will be prepared from the first element of the stride array **/
	int cur_index_init = 1;		// initialising current index for generating other elements of stride array 
	while (true)	{
		if (cur_index_init == basis_size_vect.size()-1)	{
			break;
	}
		else 
		{
			cur_product = int(cur_index_pdt/basis_size_vect[cur_index_init]);
			cur_index_init +=1;
			stride_arr_init.emplace_back(cur_product);
			cur_index_pdt = cur_product;
		}
	}
	return stride_arr_init;
}


/**********************************************************
 ********* member function to generate multidimensional 
			index from onedimensional index array  *******
 **********************************************************/

int STRIDE_ARR_BACK_N_FORTH::multidim_index_dpb(const std::vector<int> &onedim_index_vect)	{
	/** Given one dimensional indices , will return multidimensional index 
	 ** array/vector containing one dimensional indices are zero-based indices **/
	std::vector<int>stride_arr = this->stride_arr();	// calling stride array
	int multidim_basis_index = onedim_index_vect.back();
	multidim_basis_index++;		// adding one for 1-based indexing 
	for (size_t i=0; i<stride_arr.size(); i++)	{
		multidim_basis_index += stride_arr[i]*onedim_index_vect[i];
	}
	return multidim_basis_index;
}


/********************************************************************************************
 * ********* member function to one dimensional index array (returned as vector)	*********
 *	*******		from multidimensional index (1-based indexing)
 *******************************************************************************************/
std::vector<int> STRIDE_ARR_BACK_N_FORTH::onedim_indices(const int & multidim_index)	{
	/** Given multidimensional index for direct product basis, returns
	 ** one dimensional indices ( 0- based indexing ); multidim_index has 1-based indexing ***/
	std::vector<int>stride_arr = this->stride_arr();	// calling stride array 
	int multidim_index_4_caln;
	multidim_index_4_caln = multidim_index -1;
	std::vector<int>onedim_index_vect;
    /** multidim_index will change for finding each of the 1-d indices in the loop .. 
	 * here it is first initialized **/
	int cur_onedim_index;
	for (size_t i=0; i<stride_arr.size(); i++)	{
		cur_onedim_index = int(multidim_index_4_caln/stride_arr[i]);
		onedim_index_vect.emplace_back(cur_onedim_index);
		multidim_index_4_caln -= cur_onedim_index*stride_arr[i];
	}
	onedim_index_vect.emplace_back(multidim_index_4_caln);
	/** returns 1 dimensional index array .. zero based indexing **/
	return onedim_index_vect;
}


/***----------------------------------------------------------***/
/***----------------------------------------------------------***/

/** function for generating displacement 
 * grids for all coordinates of interest for gradient calculation 
 *  using 5-Pt numerical differentiation **/

std::vector<std::vector<double>> get_all_onedim_grids_4_grad(const std::vector<double>& curgrid)	{
//	std::vector<std::vector<double>> multgrid_4_numer_drv(curgrid.size(),std::vector<double>(5,0.));
	std::vector<std::vector<double>> multgrid_4_numer_drv(SD_PARAMS::mode_id_2_optimise.size(),std::vector<double>(5,0.));
	for (size_t i=0; i<multgrid_4_numer_drv.size(); i++)	{
		multgrid_4_numer_drv[i][0] = curgrid[SD_PARAMS::mode_id_2_optimise[i]-1]-2*SD_PARAMS::drv_stepsize;
		multgrid_4_numer_drv[i][1] = curgrid[SD_PARAMS::mode_id_2_optimise[i]-1]-SD_PARAMS::drv_stepsize;
		multgrid_4_numer_drv[i][2] = curgrid[SD_PARAMS::mode_id_2_optimise[i]-1];
		multgrid_4_numer_drv[i][3] = curgrid[SD_PARAMS::mode_id_2_optimise[i]-1]+SD_PARAMS::drv_stepsize;
		multgrid_4_numer_drv[i][4] = curgrid[SD_PARAMS::mode_id_2_optimise[i]-1]+2*SD_PARAMS::drv_stepsize;
	}
	return multgrid_4_numer_drv;
}




/** function to generate multidimensional grids 
 * for function evaluation at displaced grids 
 * for all coordinates to be used for gradient 
 * calculation with 5-Pt numerical differentiation ****/

std::vector<std::vector<double>> get_multidim_grid_from_all_onedim_grids(const std::vector<std::vector<double>>& all_onedim_grids,
		const std::vector<double>& const_coord_vals)	{
	std::vector<std::vector<double>> multigrid_4_grad;
	for (size_t i=0; i<all_onedim_grids.size(); i++)	{
		for (size_t j=0; j<all_onedim_grids[i].size(); j++)	{
			std::vector<double> cur_onedim_vect = const_coord_vals;
			cur_onedim_vect[SD_PARAMS::mode_id_2_optimise[i]-1] = all_onedim_grids[i][j];
			multigrid_4_grad.push_back(cur_onedim_vect);
			cur_onedim_vect.clear();
		}
	}
	return multigrid_4_grad;
}

