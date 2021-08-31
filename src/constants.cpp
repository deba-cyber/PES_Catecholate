
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
	extern const std::string gridfile = "../input_data/contracted_setB_trial_3";
	extern const std::string pot_savefile = "../output_data/Potential_setB_trial_3";
	extern const std::vector<double> alpha_vect = {6.9, 2.58, 2.60, 2.50, 2.0, 2.0};
}

namespace SD_PARAMS {
    extern const std::vector<int> mode_id_vect = {1,5,7,10,13,27,29};
    extern const std::vector<int> mode_id_2_optimise = {1,5,7,10,13,27,29};
    extern const std::vector<double> startvect = {-2.0,0.5,0.1,2.8,-0.1,0.1,1.0};
}


