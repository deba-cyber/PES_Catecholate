
# include <PES_catecholate/PES.hpp>
# include <PES_catecholate/constants.hpp>
# include <PES_catecholate/fileops.hpp>
# include <iostream>
# include <vector>
# include <string>





int main()	{
	/*** dummyvectors of int & double types  .. will be used whenever needed ***/
	std::vector<double> dummyvec_double;
	std::vector<int> dummyvec_int;
	/**** necessary parameters for PES_catecholate class ***/
	double prefact_4_scaled_normal_coords = PES_PARAMS::prefact;
	std::vector<int> Mode_flag_TS = get_specific_line(1,dummyvec_int,PES_PARAMS::modeflagfile);
	std::vector<int> linenovect_4_ref_ts = Get_lineno_vect(PES_PARAMS::N_lines_ref_ts);
	std::vector<int> linenovect_4_ref_eqm = Get_lineno_vect(PES_PARAMS::N_lines_ref_eqm);
	std::vector<int> linenovect_4_mass = Get_lineno_vect(PES_PARAMS::N_lines_mass);
	std::vector<int> linenovect_4_l_TS = Get_lineno_vect(PES_PARAMS::N_lines_l_TS);
	std::vector<int> linenovect_4_eqm_rot_mat = Get_lineno_vect(PES_PARAMS::N_lines_eqm_rot_mat);
	std::vector<int> linenovect_4_l_eqm = Get_lineno_vect(PES_PARAMS::N_lines_l_eqm);
	std::vector<int> linenovect_4_ratio_TS = Get_lineno_vect(PES_PARAMS::N_lines_ratio_ts);
	std::vector<int> linenovect_4_freq_EQM = Get_lineno_vect(PES_PARAMS::N_lines_freq_EQM);
	std::vector<double> REF_TS = get_all_sp_lines(linenovect_4_ref_ts,dummyvec_double,
			PES_PARAMS::ref_ts_file);
	std::vector<double> REF_EQM = get_all_sp_lines(linenovect_4_ref_eqm,dummyvec_double,
			PES_PARAMS::ref_eqm_file);
	std::vector<double> l_TS = get_all_sp_lines(linenovect_4_l_TS,dummyvec_double,
			PES_PARAMS::l_TS_file);
	std::vector<double> EQM_rot_mat = get_all_sp_lines(linenovect_4_eqm_rot_mat,dummyvec_double,
			PES_PARAMS::eqm_rot_mat_file);
	std::vector<double> l_eqm = get_all_sp_lines(linenovect_4_l_eqm,dummyvec_double,
			PES_PARAMS::l_eqm_file);
	std::vector<double> mass_sqrt = get_all_sp_lines(linenovect_4_mass,dummyvec_double,
			PES_PARAMS::massfile);
	std::vector<double> ratio_TS = get_all_sp_lines(linenovect_4_ratio_TS,dummyvec_double,
			PES_PARAMS::ratio_ts_file);
	std::vector<double> freq_EQM = get_all_sp_lines(linenovect_4_freq_EQM,dummyvec_double,
			PES_PARAMS::freq_eqm_file); 
	/*** Preparing Q_K_TS_scaled matrix (2-dimensional vector) ***/
	std::vector<std::vector<double>> Q_K_TS_scaled = {6,std::vector<double>(33,0.)};
	Q_K_TS_scaled[1][9] = 1.5;
	Q_K_TS_scaled[2][9] = 3.0;
	Q_K_TS_scaled[3][9] = 4.0;
	Q_K_TS_scaled[4][0] = -1.0;
	Q_K_TS_scaled[4][9] = 3.25;
	Q_K_TS_scaled[5][0] = 1.0;
	Q_K_TS_scaled[5][9] = 3.25;
	/*** Preparing Coefficient vector for coupling term ****/
	std::vector<std::vector<std::vector<double>>> B_fit_vect;
	std::vector<int> lineno_vect_4_B_vect = Get_lineno_vect(PES_PARAMS::N_lines_B_vect);
	for (int i=0; i<6; i++)	{
		std::string cur_filestr = PES_PARAMS::B_vect_file+std::to_string(i);
		std::vector<std::vector<double>> cur_B_vect_cont = Get_all_lines_multvect(lineno_vect_4_B_vect,dummyvec_double,cur_filestr);
		B_fit_vect.push_back(cur_B_vect_cont);
	}
	/*******************************************************/
	/** Creating PES_catecholate object **/
	PES_catecholate pesobj(prefact_4_scaled_normal_coords,PES_PARAMS::V_eq,PES_PARAMS::morse_a,
			Mode_flag_TS,REF_TS,REF_EQM,l_TS,EQM_rot_mat,l_eqm,mass_sqrt,ratio_TS,
			Q_K_TS_scaled,freq_EQM,B_fit_vect,PES_PARAMS::alpha_vect);
	PES_catecholate* pesobjptr = &pesobj;	// Pointer to pesobj
	/*******************************************************/
	/*******************************************************/
	/****************** Creating Diabat object ************/
	std::vector<int> linenovect_4_diabat_quad_coeff = Get_lineno_vect(PES_PARAMS::N_lines_quad_coeff); 
	std::vector<int> linenovect_4_diabat_cub_coeff = Get_lineno_vect(PES_PARAMS::N_lines_cub_coeff);
	std::vector<int> linenovect_4_diabat_quartic_coeff = Get_lineno_vect(PES_PARAMS::N_lines_quartic_coeff);
	std::vector<int> linenovect_4_diabat_high_crr = Get_lineno_vect(PES_PARAMS::N_lines_high_crr);
	std::vector<std::vector<double>> final_Morse_fact_quad = Get_all_lines_multvect(linenovect_4_diabat_quad_coeff,
			dummyvec_double,PES_PARAMS::diabat_quad_coeff_file);
	std::vector<std::vector<double>> final_Morse_fact_cub = Get_all_lines_multvect(linenovect_4_diabat_cub_coeff,
			dummyvec_double,PES_PARAMS::diabat_cub_coeff_file);
	std::vector<std::vector<double>> final_Morse_fact_quartic = Get_all_lines_multvect(linenovect_4_diabat_quartic_coeff,
			dummyvec_double,PES_PARAMS::diabat_quartic_coeff_file);
	std::vector<std::vector<double>> QS_CRR_OFF_DIAG = Get_all_lines_multvect(linenovect_4_diabat_high_crr,
			dummyvec_double,PES_PARAMS::diabat_high_crr_file);
	Diabat diabatobj(final_Morse_fact_quad,final_Morse_fact_cub,final_Morse_fact_quartic,
			QS_CRR_OFF_DIAG,PES_PARAMS::SLOPE_OO,PES_PARAMS::INTERCEPT_OO,
			PES_PARAMS::SLOPE_C2H9,PES_PARAMS::INTERCEPT_C2H9);
	/*******************************************************/
	/*******************************************************/
	std::vector<int> linenovect_4_grid = Get_lineno_vect(PES_PARAMS::N_pts);
	std::vector<std::vector<double>> multidim_grid = Get_all_lines_multvect(linenovect_4_grid,
			dummyvec_double,PES_PARAMS::gridfile); 
	std::vector<double> E_adiabat_all_vect(multidim_grid.size());
	for (size_t i=0; i<multidim_grid.size(); i++)	{
		std::vector<double> Q_TS_scaled_cur = multidim_grid[i];
		std::vector<double> V11_V22 = diabatobj.Get_tot_V11_and_V22(Q_TS_scaled_cur,pesobjptr);
		double V12_sqr = pesobj.V12_sqr(Q_TS_scaled_cur);
		double E_adiabat_cur = Get_Adiabat(V11_V22,V12_sqr);
		E_adiabat_all_vect[i] = E_adiabat_cur;
	}
	save_binary(E_adiabat_all_vect,PES_PARAMS::pot_savefile,"append");
	return EXIT_SUCCESS;
}

