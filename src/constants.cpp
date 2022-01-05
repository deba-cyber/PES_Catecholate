
# include <PES_catecholate/constants.hpp>
# include <iosfwd>
# include <vector>
# include <string>

namespace PES_PARAMS	{
	extern const std::string modeflagfile = "../input_data/mode_flag_TS";
	extern const std::string ref_ts_file = "../input_data/ref_com_rmv_TS";
	extern const std::string ref_eqm_file = "../input_data/ref_com_rmv_eqm";
	extern const std::string l_TS_file = "../input_data/catecholate_TS_l-matrix";
	extern const std::string eqm_rot_mat_file = "../input_data/Eckart_mat_eqm";
	extern const std::string l_eqm_file = "../input_data/l_eqm_muw_rnrm";	
	extern const std::string massfile = "../input_data/mass_sqrt";
	extern const std::string ratio_ts_file = "../input_data/ratio_TS";
	extern const std::string B_vect_file = "../input_data/Q_fitparams_";
	extern const std::string freq_eqm_file = "../input_data/freq_eqm";
	extern const std::string diabat_quad_coeff_file = "../input_data/chk_morse_hessian_only_diag_1";
	extern const std::string diabat_cub_coeff_file = "../input_data/chk_morse_com_ref_cubic_9_pt";
	extern const std::string diabat_quartic_coeff_file = "../input_data/chk_morse_com_ref_quartic_9_pt";
	extern const std::string diabat_high_crr_file = "../input_data/QS_CRR_OFF_DIAG";
	extern const std::string gridfile = "../input_data/E_IM_grid_RMS";
	extern const std::string pot_savefile = "../output_data/Potential_IM_grid";
	extern const std::vector<double> alpha_vect = {6.9, 2.58, 2.60, 2.50, 2.0, 2.0};
	extern const std::string datasavetype = "dat";
}

namespace SD_PARAMS	{
	extern const std::string gridfile_2_optimise = "../input_data/Q1_Q10_relax_input_grid";
	extern const std::string gridfile_2_start = "../input_data/Q1_Q10_relax_starvect";
	extern const std::vector<int> mode_id_vect = {1,5,7,8,10,13,16,17,18,19,21,22,23,24,25,26,27,28,29,30,31,32,33};
	extern const std::vector<int> mode_id_2_optimise = {1,5,7,8,13,16,17,18,19,21,22,23,24,25,26,27,28,30,31,32,33};
//	extern const std::vector<int> mode_id_fixed = {1,5,10,29};
//	extern const std::vector<double> startvect = {-2.0,0.5,0.3,2.8,-0.3,-0.116141,0.202347,0.152987,0.0969398,-0.225726,0.1,1.0};
	extern const std::vector<double> startvect = {-2.86,0.55,0.34,0.07,2.77,-0.34,0.12,-0.09,-0.09,0.0075,-0.1942,-0.1485,0.1013,-0.2284,0.0088,-0.0952,0.2884,-0.0361,0.83,0.0034,-0.0003,0.0052,0.001};
//	extern const std::vector<double> startvect = {7.180125e-01,5.605201e-01,-1.877045e-01,-5.937518e-01,
//		-1.646329e-01,-7.822349e-02,9.285675e-02,-4.417326e-02,3.478453e-01,1.164222e-01,2.681571e-01,-3.813042e-01,
//		5.104468e-02,-1.897087e-01,4.958839e-01,-6.497944e-02,1.457136e+00,-2.610049e-03,-8.922957e-03,9.870465e-04,
//		-4.750531e-03};	
}
