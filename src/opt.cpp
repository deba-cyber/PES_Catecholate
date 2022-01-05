# include <PES_catecholate/PES.hpp>
# include <PES_catecholate/constants.hpp>
# include <PES_catecholate/fileops.hpp>
# include <PES_catecholate/SD.hpp>
# include <iostream>
# include <cmath>
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
	/**************** Optimisation ************************/
	std::vector<double> X0_fulldim(PES_PARAMS::N_vibs);		// location to be updated in each iteration
	for (size_t i=0; i<SD_PARAMS::mode_id_vect.size(); i++)	{
		X0_fulldim[SD_PARAMS::mode_id_vect[i]-1] = SD_PARAMS::startvect[i];
	}
	std::vector<double> cur_grad_vect;
	std::vector<double> cur_grad_vect_fulldim(PES_PARAMS::N_vibs);
	{
		std::vector<std::vector<double>> all_onedim_dsplcmnts = get_all_onedim_grids_4_grad(X0_fulldim); 
		std::vector<std::vector<double>> multgrid_4_grad = get_multidim_grid_from_all_onedim_grids(all_onedim_dsplcmnts,X0_fulldim);
		std::vector<double> multidim_grid_4_func_eval;
		std::vector<std::vector<double>> func_4_grad;
		for (size_t i=0; i<multgrid_4_grad.size(); i++)	{
			std::vector<double> Q_TS_scaled_cur = multgrid_4_grad[i];
			std::vector<double> V11_V22 = diabatobj.Get_tot_V11_and_V22(Q_TS_scaled_cur,pesobjptr);
			double V12_sqr = pesobj.V12_sqr(Q_TS_scaled_cur);
			double E_adiabat_cur = Get_Adiabat(V11_V22,V12_sqr);
			multidim_grid_4_func_eval.emplace_back(E_adiabat_cur);	
		}
		func_4_grad = split_onedim_vector(multidim_grid_4_func_eval,5);
		cur_grad_vect = get_first_drv(func_4_grad,SD_PARAMS::drv_stepsize);
		for (size_t i=0; i<cur_grad_vect.size(); i++)	{
			cur_grad_vect_fulldim[SD_PARAMS::mode_id_2_optimise[i]-1] = cur_grad_vect[i];
		}
	}
	double grad_norm = Vec_Norm(cur_grad_vect_fulldim);
	std::vector<double> cur_direction_vect;
	int counter = 0;
	double alpha_cur_linsrch;
	std::vector<double> cur_test_loc_4_alpha;
	std::vector<double> tmp_direction_vect_4_alpha;
	while (grad_norm > SD_PARAMS::grad_thr && counter < SD_PARAMS::MAXITR && alpha_cur_linsrch > SD_PARAMS::alpha_thr)	{
		int alpha_ctr = 1;
		alpha_cur_linsrch = SD_PARAMS::alpha_init_4_linsrch;
		while (alpha_ctr != 0)	{
			cur_test_loc_4_alpha = X0_fulldim;
			tmp_direction_vect_4_alpha = cur_grad_vect_fulldim;
			multiply_w_const(tmp_direction_vect_4_alpha,-1.0/grad_norm);
			multiply_w_const(tmp_direction_vect_4_alpha,alpha_cur_linsrch);
			add_2_vect(cur_test_loc_4_alpha,tmp_direction_vect_4_alpha);
			double E_adiabat_cur_4_alpha;
			double E_adiabat_prv_4_alpha;
			{
				std::vector<double> V11_V22 = diabatobj.Get_tot_V11_and_V22(cur_test_loc_4_alpha,pesobjptr);
				std::vector<double> V11_V22_prv = diabatobj.Get_tot_V11_and_V22(X0_fulldim,pesobjptr);
				double V12_sqr = pesobj.V12_sqr(cur_test_loc_4_alpha);
				double V12_sqr_prv = pesobj.V12_sqr(X0_fulldim); 
				E_adiabat_cur_4_alpha = Get_Adiabat(V11_V22,V12_sqr);	
				E_adiabat_prv_4_alpha = Get_Adiabat(V11_V22_prv,V12_sqr_prv);
			}
			if (E_adiabat_cur_4_alpha > E_adiabat_prv_4_alpha+((pow(grad_norm,2.))*alpha_cur_linsrch*SD_PARAMS::beta_4_linsrch))	{
				alpha_cur_linsrch *= SD_PARAMS::tau_4_linsrch;
			}
			else
				alpha_ctr = 0;
		}
		cur_direction_vect = cur_grad_vect_fulldim;
		multiply_w_const(cur_direction_vect,-1.0/grad_norm);
		multiply_w_const(cur_direction_vect,alpha_cur_linsrch);
		add_2_vect(X0_fulldim,cur_direction_vect);
		std::vector<std::vector<double>> all_onedim_dsplcmnts = get_all_onedim_grids_4_grad(X0_fulldim);	
		std::vector<std::vector<double>> multgrid_4_grad = get_multidim_grid_from_all_onedim_grids(all_onedim_dsplcmnts,X0_fulldim);
		std::vector<double> multidim_grid_4_func_eval;
		std::vector<std::vector<double>> func_4_grad;
		for (size_t i=0; i<multgrid_4_grad.size(); i++)	{
			std::vector<double> Q_TS_scaled_cur = multgrid_4_grad[i]; 
			std::vector<double> V11_V22 = diabatobj.Get_tot_V11_and_V22(Q_TS_scaled_cur,pesobjptr);
			double V12_sqr = pesobj.V12_sqr(Q_TS_scaled_cur);
			double E_adiabat_cur = Get_Adiabat(V11_V22,V12_sqr);
			multidim_grid_4_func_eval.emplace_back(E_adiabat_cur);
		}
		func_4_grad = split_onedim_vector(multidim_grid_4_func_eval,5);
		cur_grad_vect = get_first_drv(func_4_grad,SD_PARAMS::drv_stepsize);
		make_vect_zero(cur_grad_vect_fulldim);
		for (size_t i=0; i<cur_grad_vect.size(); i++)	{
			cur_grad_vect_fulldim[SD_PARAMS::mode_id_2_optimise[i]-1] = cur_grad_vect[i];
		}
		grad_norm = Vec_Norm(cur_grad_vect_fulldim);
		counter++;
		std::cout	<<	" counter "	<< counter	<< " alpha "	<< alpha_cur_linsrch	<< " grad norm	"	<< grad_norm	<< "\n"; 
	}
	std::vector<double> V11_V22_opt = diabatobj.Get_tot_V11_and_V22(X0_fulldim,pesobjptr);
	double V12_sqr_opt = pesobj.V12_sqr(X0_fulldim);
	double E_adiabat_opt = Get_Adiabat(V11_V22_opt,V12_sqr_opt);
	std::cout	<< X0_fulldim[0]	<< "	"	<< X0_fulldim[4]	<< "	"	<< X0_fulldim[6]	<<	
		"	"	<< X0_fulldim[9]	<<	"	"	<< X0_fulldim[12]	<< "	"	<< X0_fulldim[26]	<< "	"
		<< X0_fulldim[28]	<< "\n";
	printonedimvect(cur_grad_vect);
	std::cout	<< E_adiabat_opt	<< "\n";
	return EXIT_SUCCESS;
}
