# include <PES_catecholate/SD.hpp>
# include <PES_catecholate/constants.hpp>
# include <iosfwd>
# include <vector>





/** function for generating displacement 
 * grids for all coordinates for gradient calculation 
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


