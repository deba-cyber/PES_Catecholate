/* file for preparing grid for reduced dimensional Hamiltonian 
 * solution in the contracted basis approach */

# include <iostream>
# include <fstream>
# include <sstream>
# include <vector>


/**************************************************
 * Function for extracting specific lines
 * from a file	.................................
 * string s1 is passed by value filename wo ext.
 * vector lineno_vect contains numbers in increasing
 * order corrsponsing to lines to be extracted	
 * lineno_vect has 1-based indexing	...............
 * returning as one dimensional vector ..........
 * specific to data type .. here double
 * ************************************************/

std::vector<double> get_all_sp_lines(const std::vector<int> & lineno_vect, const std::string& s1)		{
	std::string line1;
	std::string ext = ".dat";
	std::string path = s1 + ext;
	std::vector<double> rtrn_multiple_lines_vect;
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
				double val; 
				while (lstream >> val)	{
					rtrn_multiple_lines_vect.emplace_back(val);
				}
			}
		}
	}
	return rtrn_multiple_lines_vect;
}


/***************************************************
 * Function for generating line number vector from 
 * number of lines ..
 * *************************************************/

std::vector<int> Get_lineno_vect(int NUMLINES)	{
	std::vector<int> linenovect(NUMLINES);
	for (size_t i=0; i<linenovect.size(); i++)	{
		linenovect[i] = i+1;
	}
	return linenovect;
}



int main()	{
	int N_po_pts_Q1 = 18;
	int N_po_pts_Q5 = 12;
	int N_po_pts_Q7 = 8;
	int N_po_pts_Q10 = 17;
	int N_po_pts_Q13 = 7; 
	int N_po_pts_Q27 = 6;
	int N_po_pts_Q29 = 7;
	std::vector<int> linenovect_Q1 = Get_lineno_vect(N_po_pts_Q1); 	
	std::vector<int> linenovect_Q5 = Get_lineno_vect(N_po_pts_Q5);
	std::vector<int> linenovect_Q7 = Get_lineno_vect(N_po_pts_Q7);
	std::vector<int> linenovect_Q10 = Get_lineno_vect(N_po_pts_Q10);
	std::vector<int> linenovect_Q13 = Get_lineno_vect(N_po_pts_Q13);
	std::vector<int> linenovect_Q27 = Get_lineno_vect(N_po_pts_Q27);
	std::vector<int> linenovect_Q29 = Get_lineno_vect(N_po_pts_Q29);
	/** 1d grids **/
	std::vector<double> Q1_grid = get_all_sp_lines(linenovect_Q1,"podvr_points_1");
	std::vector<double> Q5_grid = get_all_sp_lines(linenovect_Q5,"podvr_points_5");
	std::vector<double> Q7_grid = get_all_sp_lines(linenovect_Q7,"podvr_points_7");
	std::vector<double> Q10_grid = get_all_sp_lines(linenovect_Q10,"podvr_points_10");
	std::vector<double> Q13_grid = get_all_sp_lines(linenovect_Q13,"podvr_points_13");
	std::vector<double> Q27_grid = get_all_sp_lines(linenovect_Q27,"podvr_points_27");
	std::vector<double> Q29_grid = get_all_sp_lines(linenovect_Q29,"podvr_points_29");
	/**--------------**/
	unsigned int DPG_3d_contracted = 8*7*6;
	unsigned int DPG_4d_contracted = 18*12*17*7;
	std::vector<std::vector<double>> grid_4_contracted_3d = {DPG_3d_contracted,std::vector<double>(33,0.)};		// 3d grid (Q7-Q13-Q27) for contracted basis calculation
	std::vector<std::vector<double>> grid_4_contracted_4d = {DPG_4d_contracted,std::vector<double>(33,0.)};		// 4d grid (Q1-Q5-Q10-Q29) for contracted basis calculation
	/**-- 3d direct product grid --**/
	int counter = 0;
	for (size_t i=0; i<Q7_grid.size(); i++)	{
		for (size_t j=0; j<Q13_grid.size(); j++)	{
			for (size_t k=0; k<Q27_grid.size(); k++)	{
			//	grid_4_contracted_3d[counter][0] = -2.0;
			//	grid_4_contracted_3d[counter][4] = -0.5;
				grid_4_contracted_3d[counter][6] = Q7_grid[i];
			//	grid_4_contracted_3d[counter][9] = 2.0;
				grid_4_contracted_3d[counter][12] = Q13_grid[j];
				grid_4_contracted_3d[counter][26] = Q27_grid[k];
			//	grid_4_contracted_3d[counter][28] = 0.8;
				counter++;
			}
		}
	}
	/**-- 4d direct product grid --**/
	int counter_4d = 0;
	for (size_t i=0; i<Q1_grid.size(); i++)	{
		for (size_t j=0; j<Q5_grid.size(); j++)	{
			for (size_t k=0; k<Q10_grid.size(); k++)	{
				for (size_t l=0; l<Q29_grid.size(); l++)	{
					grid_4_contracted_4d[counter_4d][0] = Q1_grid[i];
					grid_4_contracted_4d[counter_4d][4] = Q5_grid[j];
					grid_4_contracted_4d[counter_4d][6] = 0.3;
					grid_4_contracted_4d[counter_4d][9] = Q10_grid[k];
					grid_4_contracted_4d[counter_4d][12] = -0.3;
					grid_4_contracted_4d[counter_4d][26] = 0.3;
					grid_4_contracted_4d[counter_4d][28] = Q29_grid[l];
					counter_4d++;
				}
			}
		}
	}
	/***----------------------------***/
	for (size_t i=0; i<grid_4_contracted_3d.size(); i++)	{
		for (size_t j=0; j<grid_4_contracted_3d[i].size(); j++)	{
			std::cout	<< std::scientific	<< grid_4_contracted_3d[i][j]	<< "	";
		}
		std::cout	<< std::endl;
	}
	return 0;
}
