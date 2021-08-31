# if ! defined (_CONST_H)
# define _CONST_H

# include <iosfwd>
# include <vector>
# include <string>


namespace PES_PARAMS	{
	inline constexpr double prefact = 0.17222125670294633;
	inline constexpr double V_eq = -1581.800133;
	inline constexpr double morse_a = 0.26672731621158;
	inline constexpr int N_vibs = 33;	// No of normal modes in the molecule
	inline constexpr int N_atoms = 13;	// No of atoms in the molecule
	/* filenames containing necessary data and constants 
	 * denoting number of lines in thise files for more than 1 line */
	extern const std::string modeflagfile;	// file containing mode flags based on symmetry
	extern const std::string ref_ts_file;	// file containing REF_TS cartesians
	inline constexpr int N_lines_ref_ts = 13;	// No of lines in ref_ts file
	extern const std::string ref_eqm_file;	// file containing REF_EQM cartesians
	inline constexpr int N_lines_ref_eqm = 13;	// No of lines in ref_eqm file
	extern const std::string l_TS_file;		// file containing TS l-matrix
	inline constexpr int N_lines_l_TS = 39;		// No of lines in TS l-matrix file	
	extern const std::string eqm_rot_mat_file;	// file containing EQM rotation matrix (for transforming  to EQM frame)
	inline constexpr int N_lines_eqm_rot_mat = 3;	// No of lines in Eqm Rotation matrix
	extern const std::string l_eqm_file;		// file containing eqm l-matrix (mass-unweighted and re-normalized)
	inline constexpr int N_lines_l_eqm = 39;	// No of lines in l_eqm file 
	extern const std::string massfile;		// file containing diagonal elements of mass^(1/2) matrix(diagonal)
	inline constexpr int N_lines_mass = 39;	// No of lines in mass file
	extern const std::string ratio_ts_file;	// file containing ratio at TS
	inline constexpr int N_lines_ratio_ts = 33;		// No of lines in ratio_ts file
	extern const std::string B_vect_file;			// string for B-vector in coupling term
	inline constexpr int N_lines_B_vect = 595;		// No of lines in Coupling coefficient vector
	extern const std::string freq_eqm_file;	// file containing eqm frequencies
	inline constexpr int N_lines_freq_EQM = 33;
	inline constexpr double SLOPE_OO = -3.272235531557268;
	inline constexpr double INTERCEPT_OO = -6.75;
	inline constexpr double SLOPE_C2H9 = -0.06579365079365072;
	inline constexpr double INTERCEPT_C2H9 = 1.128;
	inline constexpr int N_lines_quad_coeff = 33;
	inline constexpr int N_lines_cub_coeff = 3565;
	inline constexpr int N_lines_quartic_coeff = 5264;
	inline constexpr int N_lines_high_crr = 9;
	extern const std::string diabat_quad_coeff_file;		// filename containing quadratic coefficients for Diabat
	extern const std::string diabat_cub_coeff_file;			// filename containing cubic coefficients for Diabat
	extern const std::string diabat_quartic_coeff_file;		// filename containinf quartic coefficients for Diabat
	extern const std::string diabat_high_crr_file;
	inline constexpr int N_pts = 336;					// Total number of multidimensional grid points 
	extern const std::string gridfile;
	extern const std::vector<double> alpha_vect;	// vector containing coefficients in the DGEVB coupling term
	extern const std::string pot_savefile;	
}

namespace SD_PARAMS	{
	inline constexpr int MAXITR = 1000;
	inline constexpr double grad_thr = 1e-02;
	inline constexpr double alpha_init_4_linsrch = 0.005;	// initial alpha for line search in each iteration
	inline constexpr double beta_4_linsrch = 0.01;			
	inline constexpr double tau_4_linsrch = 0.5;
	inline constexpr double first_drv_coeff_1 = 0.833333333333333E-01; 
	inline constexpr double first_drv_coeff_2 = -0.666666666666667E+00;
	inline constexpr double first_drv_coeff_3 = 0.594762334620620E-16;
	inline constexpr double first_drv_coeff_4 = 0.666666666666667E+00;
	inline constexpr double first_drv_coeff_5 = -0.833333333333333E-01;
	inline constexpr double drv_stepsize = 0.01;
	extern const std::vector<int> mode_id_vect;
	extern const std::vector<int> mode_id_2_optimise;
	extern const std::vector<double> startvect;
}

# endif
