
# include <PES_catecholate/PES.hpp>
# include <PES_catecholate/fileops.hpp>
# include <PES_catecholate/constants.hpp>
# include <iostream>
# include <vector>
# include <cmath>

/****************************************
 * Class constructor of PES_catecholate
 * **************************************/


PES_catecholate::PES_catecholate(double prefact_4_scaled_normal_coords,double V_eq,double morse_a,
		std::vector<int>& Mode_flag_TS,std::vector<double>& REF_TS,
		std::vector<double>& REF_EQM,std::vector<double>& l_TS,
		std::vector<double>& EQM_rot_mat,std::vector<double>& l_eqm,
		std::vector<double>& mass_sqrt,std::vector<double>& ratio_TS,
		std::vector<std::vector<double>>& Q_K_TS_scaled,std::vector<double>& freq_EQM,
		std::vector<std::vector<std::vector<double>>>& B_fit_vect,
		const std::vector<double>& alpha_vect)	{
	this-> prefact_4_scaled_normal_coords = prefact_4_scaled_normal_coords;
	this-> V_eq = V_eq;
	this-> morse_a = morse_a;
	this-> Mode_flag_TS = Mode_flag_TS;
	this-> REF_EQM = REF_EQM;
	this-> REF_TS = REF_TS;
	this-> l_TS = l_TS;
	this-> EQM_rot_mat = EQM_rot_mat;
	this-> l_eqm = l_eqm;
	this-> mass_sqrt = mass_sqrt;
	this-> ratio_TS = ratio_TS;
	this-> Q_K_TS_scaled = Q_K_TS_scaled;
	this-> freq_EQM = freq_EQM;
	this-> B_fit_vect = B_fit_vect;
	this-> alpha_vect = alpha_vect;
}


/*************************************************************************
 * Member function to transform scaled coordinates(TS) to Gauss units
 * **********************************************************************/


std::vector<double> PES_catecholate::TS_scaled_2_Gauss(const std::vector<double>& Q_TS_scl) const {
	std::vector<double> Q_TS_Gauss(Q_TS_scl.size());
	for (size_t i=0; i<Q_TS_Gauss.size(); i++)	{
		if (this->ratio_TS[i] != 0.)	{
			Q_TS_Gauss[i] = Q_TS_scl[i]/this->ratio_TS[i];
		}
	}
	return Q_TS_Gauss;
}


/*************************************************************************
 * Member function to transform EQM renormalized normal coordinates to 
 * corresponding dimension-less normal coordinates 
 * ***********************************************************************/


std::vector<double> PES_catecholate::Eqm_rnrm_2_scaled(const std::vector<double>& Q_eq_rnrm) const {
	std::vector<double> Q_eq_scaled(Q_eq_rnrm.size());
	for (size_t i=0; i<Q_eq_scaled.size(); i++)	{
		Q_eq_scaled[i] = (this->prefact_4_scaled_normal_coords)*
			sqrt(std::abs(this->freq_EQM[i]))*Q_eq_rnrm[i];		
	}
	return Q_eq_scaled;
}


/*************************************************************************
 * Member function to calculate V12 square corresponding to a multidim
 * grid in TS normal coordinates
*************************************************************************/

double PES_catecholate::V12_sqr(const std::vector<double>& Q_TS_scaled_cur) const {
	double V12_sqr_tot = 0.;
	for (size_t i=0; i<this->B_fit_vect.size(); i++)	{
		std::vector<double> Q_diff_vect = subtract_vect(Q_TS_scaled_cur,this->Q_K_TS_scaled[i]); 
		double norm_sqr = get_norm_sqr(Q_diff_vect);
		double G0 = exp(-0.5*(this->alpha_vect[i])*norm_sqr);
		for (size_t j=0; j<this->B_fit_vect[i].size(); j++)	{
			if (this->B_fit_vect[i][j][0] == this->B_fit_vect[i][j][1] && static_cast<int>(this->B_fit_vect[i][j][1]) == 0)	{
				V12_sqr_tot += this->B_fit_vect[i][j][2]*G0; 
			}
			else if (this->B_fit_vect[i][j][1] == 0)	{
				V12_sqr_tot += this->B_fit_vect[i][j][2]*(Q_TS_scaled_cur[static_cast<int>(this->B_fit_vect[i][j][0])-1]-Q_K_TS_scaled[i][static_cast<int>(this->B_fit_vect[i][j][0])-1])
					*G0;
			}
			else if (this->B_fit_vect[i][j][0] == 0)	{
				V12_sqr_tot += this->B_fit_vect[i][j][2]*(Q_TS_scaled_cur[static_cast<int>(this->B_fit_vect[i][j][1])-1]-Q_K_TS_scaled[i][static_cast<int>(this->B_fit_vect[i][j][1])-1])
					*G0;
			}
			else {
				V12_sqr_tot += this->B_fit_vect[i][j][2]*(Q_TS_scaled_cur[static_cast<int>(this->B_fit_vect[i][j][0])-1]-Q_K_TS_scaled[i][static_cast<int>(this->B_fit_vect[i][j][0])-1])
					*(Q_TS_scaled_cur[static_cast<int>(this->B_fit_vect[i][j][1])-1]-Q_K_TS_scaled[i][static_cast<int>(this->B_fit_vect[i][j][1])-1])*G0;
			}
		}
	}
	return V12_sqr_tot;
}


/*************************************
 * Class constructor to class Diabat
**************************************/


Diabat::Diabat(std::vector<std::vector<double>>& final_Morse_fact_quad,std::vector<std::vector<double>>& final_Morse_fact_cub,
		std::vector<std::vector<double>>& final_Morse_fact_quartic,std::vector<std::vector<double>>& QS_CRR_OFF_DIAG,
		double mean_slope_OO,double intercept_OO,double mean_slope_C2H9,double intercept_C2H9)	{
			this-> final_Morse_fact_quad = final_Morse_fact_quad;
			this-> final_Morse_fact_cub = final_Morse_fact_cub;
			this-> final_Morse_fact_quartic = final_Morse_fact_quartic;
			this-> QS_CRR_OFF_DIAG = QS_CRR_OFF_DIAG;
			this-> mean_slope_OO = mean_slope_OO;
			this-> intercept_OO = intercept_OO;
			this-> mean_slope_C2H9 = mean_slope_C2H9;
			this-> intercept_C2H9 = intercept_C2H9;
}

/*******************************************************************
 * Member function to calculate Mirror grid corresponding to given
 * TS multidimensional grid ......................................
********************************************************************/

std::vector<double> Diabat::Get_Q_TS_mirror_grid(const std::vector<double>& Q_TS_scaled_cur,PES_catecholate* pesobj) const {
	std::vector<double> Q_TS_mirror(Q_TS_scaled_cur.size());
	for (size_t i=0; i<Q_TS_mirror.size(); i++)	{
		if (pesobj->Mode_flag_TS[i] !=	1)	{
			Q_TS_mirror[i] = -Q_TS_scaled_cur[i];
		}
		else	{
			Q_TS_mirror[i] = Q_TS_scaled_cur[i];
		}
	}
	return Q_TS_mirror;
}

/******************************************************************
 * Member function to generate Q_eqm normal coordinates from Q_TS 
 * coordinates .. If Q_TS_mirror grid is passed, it will generate
 * Q_eqm_mirror grid 
 * ****************************************************************/


std::vector<double> Diabat::Get_Q_EQM_4_TS_or_Mirror(const std::vector<double>& Q_TS_scaled_cur,
		PES_catecholate* pesobj) const	{
	std::vector<double> Q_TS_Gauss = pesobj->TS_scaled_2_Gauss(Q_TS_scaled_cur);
	std::vector<double> delta_cart_TS = get_mat_vec_pdt(pesobj->l_TS,Q_TS_Gauss,3*(PES_PARAMS::N_atoms),PES_PARAMS::N_vibs);
	std::vector<double> cart_config_TS(delta_cart_TS.size());
	for (size_t i=0; i<cart_config_TS.size(); i++)	{
		cart_config_TS[i] = delta_cart_TS[i] + pesobj->REF_TS[i];
	}
	// cart_config_TS gives cartesians in TS frame //
	// ---------------------------------------------//
	// calculating Q_eqm(normal modes at eqm w.r.t. EQM frame) //
	std::vector<double> cart_config_EQM;	// cart_config_EQM gives cartesians in EQM frame
	for (size_t i=0; i<3*(PES_PARAMS::N_atoms); i+=3)	{
		std::vector<double> cur_cart = {cart_config_TS[i],cart_config_TS[i+1],cart_config_TS[i+2]};
		std::vector<double> tmp_rotated_vect = get_mat_vec_pdt(pesobj->EQM_rot_mat,cur_cart,3,3);
		cart_config_EQM.insert(std::end(cart_config_EQM),std::begin(tmp_rotated_vect),std::end(tmp_rotated_vect));
		empty_swap(cur_cart);
		empty_swap(tmp_rotated_vect);
	}
	// cartesian displacement vector wrt EQM frame //
	std::vector<double> delta_cart_EQM = subtract_vect(cart_config_EQM,pesobj->REF_EQM);
	// multiplying delta_cart_EQM with m^(1/2) matrix //
	for (size_t i=0; i<delta_cart_EQM.size(); i++)	{
		delta_cart_EQM[i] *= pesobj->mass_sqrt[i];
	}
	// getting transpose of mass-unweighted renormalized EQM l-matrix as 1-dim vector //
	std::vector<double> l_eqm_T = get_mat_transpose_as_vect(pesobj-> l_eqm,3*(PES_PARAMS::N_atoms),PES_PARAMS::N_vibs); 
	// getting Q_eqm_rnrm //
	std::vector<double> Q_eqm_rnrm = get_mat_vec_pdt(l_eqm_T,delta_cart_EQM,PES_PARAMS::N_vibs,3*(PES_PARAMS::N_atoms)); 
	// getting Q_scaled eqm //
	std::vector<double> Q_eqm_scaled = pesobj->Eqm_rnrm_2_scaled(Q_eqm_rnrm);
	for (size_t i=0; i<Q_eqm_scaled.size(); i++)	{
		if (std::abs(Q_eqm_scaled[i]) < 5e-04)	{
			Q_eqm_scaled[i] = 0.0;
		}
	}
	return Q_eqm_scaled;
}


/******************************************************************
 * Member function to calculate Morse coordinate corresponding 
 * to multidimensional Eqm. grid .............................
******************************************************************/


std::vector<double> Diabat::Get_Q_eqm_morse_grid(const std::vector<double>& Q_eqm_grid,PES_catecholate* pesobj) const {
	std::vector<double> Q_eqm_morse_grid(Q_eqm_grid.size());
	for (size_t i=0; i<Q_eqm_morse_grid.size()-1; i++)	{
		Q_eqm_morse_grid[i] = Q_eqm_grid[i];
	}
	Q_eqm_morse_grid[Q_eqm_morse_grid.size()-1] = 1-exp(-pesobj->morse_a*Q_eqm_grid[Q_eqm_grid.size()-1]);
	return Q_eqm_morse_grid;
}

/****************************************************************************
 * Member function for returning relevant internal parameters and EQM
 * normal mode displacements corresponding to a particular TS configuration
 * Returning a 2-dimensional vector .... rtrn_internal_vect
 * *************************************************************************/

std::vector<std::vector<double>> Diabat::Get_all_internal_params(const std::vector<double>& Q_TS_scaled_cur,
		PES_catecholate* pesobj) const	{
	std::vector<std::vector<double>> rtrn_internal_vect = {{0.},{0.,0.},
		{0.,0.},{0.,0.,0.,0.},{0.}};
	std::vector<double> Q_TS_Gauss = pesobj->TS_scaled_2_Gauss(Q_TS_scaled_cur);
	std::vector<double> delta_cart_TS = get_mat_vec_pdt(pesobj->l_TS,Q_TS_Gauss,3*(PES_PARAMS::N_atoms),PES_PARAMS::N_vibs);
	std::vector<double> cart_config_TS(delta_cart_TS.size());
	// cart_config_TS gives cartesians in TS frame //
	for (size_t i=0; i<cart_config_TS.size(); i++)	{
		cart_config_TS[i] = delta_cart_TS[i] + pesobj->REF_TS[i];
	}
	/******************************************************/
	/* calculating r_OH  i.e. distance between O4 and H9 */
	rtrn_internal_vect[0][0] = sqrt(pow((cart_config_TS[25]-cart_config_TS[10]),2.)+pow((cart_config_TS[26]-cart_config_TS[11]),2.));
	/* calculating O-ortho-H repulsion parameters .. r_O1H13,r_O4H10 */
	rtrn_internal_vect[1][0] = sqrt(pow((cart_config_TS[37]-cart_config_TS[1]),2.)+pow((cart_config_TS[38]-cart_config_TS[2]),2.));
	rtrn_internal_vect[1][1] = sqrt(pow((cart_config_TS[28]-cart_config_TS[10]),2.)+pow((cart_config_TS[29]-cart_config_TS[11]),2.));
	/******************************************************/
	/* calculating O-O repulsion parameters .. aO1C2C3,aO4C3C2 */
	double r_C2C3 = sqrt(pow((cart_config_TS[7]-cart_config_TS[4]),2.)+pow((cart_config_TS[8]-cart_config_TS[5]),2.));
	double r_C2O1 = sqrt(pow((cart_config_TS[1]-cart_config_TS[4]),2.)+pow((cart_config_TS[2]-cart_config_TS[5]),2.));
	double r_C3O4 = sqrt(pow((cart_config_TS[10]-cart_config_TS[7]),2.)+pow((cart_config_TS[11]-cart_config_TS[8]),2.));
	std::vector<double> vec_C3O4 = {cart_config_TS[10]-cart_config_TS[7],cart_config_TS[11]-cart_config_TS[8]};
	std::vector<double> vec_C3C2 = {cart_config_TS[4]-cart_config_TS[7],cart_config_TS[5]-cart_config_TS[8]};
	std::vector<double> vec_C2O1 = {cart_config_TS[1]-cart_config_TS[4],cart_config_TS[2]-cart_config_TS[5]};
	std::vector<double> vec_C2C3 = {cart_config_TS[7]-cart_config_TS[4],cart_config_TS[8]-cart_config_TS[5]};
	double aO1C2C3 = std::acos((get_dot_pdt(vec_C2O1,vec_C2C3))/(r_C2C3*r_C2O1))*180.*(1/M_PI);
	double aO4C3C2 = std::acos((get_dot_pdt(vec_C3O4,vec_C3C2))/(r_C3O4*r_C2C3))*180.*(1/M_PI);
	rtrn_internal_vect[2][0] = aO1C2C3;
	rtrn_internal_vect[2][1] = aO4C3C2;
	/******************************************************/
	/* calculating H9-C2C3 repulsion parameters .. r_C2H9,r_C3H9,aC2O1H9,aC3O4H9 */
	rtrn_internal_vect[3][0] = sqrt(pow((cart_config_TS[25]-cart_config_TS[4]),2.)+pow((cart_config_TS[26]-cart_config_TS[5]),2.));
	rtrn_internal_vect[3][1] = sqrt(pow((cart_config_TS[25]-cart_config_TS[7]),2.)+pow((cart_config_TS[26]-cart_config_TS[8]),2.));
	double r_O1C2 = sqrt(pow((cart_config_TS[4]-cart_config_TS[1]),2.)+pow((cart_config_TS[5]-cart_config_TS[2]),2.));
	double r_O1H9 = sqrt(pow((cart_config_TS[25]-cart_config_TS[1]),2.)+pow((cart_config_TS[26]-cart_config_TS[2]),2.));
	double r_O4C3 = sqrt(pow((cart_config_TS[7]-cart_config_TS[10]),2.)+pow((cart_config_TS[8]-cart_config_TS[11]),2.));
	double r_O4H9 = sqrt(pow((cart_config_TS[25]-cart_config_TS[10]),2.)+pow((cart_config_TS[26]-cart_config_TS[11]),2.));
	std::vector<double> vec_O1C2 = {cart_config_TS[4]-cart_config_TS[1],cart_config_TS[5]-cart_config_TS[2]};
	std::vector<double> vec_O1H9 = {cart_config_TS[25]-cart_config_TS[1],cart_config_TS[26]-cart_config_TS[2]};
	std::vector<double> vec_O4C3 = {cart_config_TS[7]-cart_config_TS[10],cart_config_TS[8]-cart_config_TS[11]};
	std::vector<double> vec_O4H9 = {cart_config_TS[25]-cart_config_TS[10],cart_config_TS[26]-cart_config_TS[11]};
	double aC2O1H9 = std::acos((get_dot_pdt(vec_O1C2,vec_O1H9))/(r_O1C2*r_O1H9))*180.*(1/M_PI);
	double aC3O4H9 = std::acos((get_dot_pdt(vec_O4C3,vec_O4H9))/(r_O4C3*r_O4H9))*180.*(1/M_PI);;
	rtrn_internal_vect[3][2] = aC2O1H9; 
	rtrn_internal_vect[3][3] = aC3O4H9;
	/******************************************************/
	/* calculating mirror r_OH parameters i.e. distance between O1 & H9 */
	rtrn_internal_vect[4][0] = sqrt(pow((cart_config_TS[25]-cart_config_TS[1]),2.)+pow((cart_config_TS[26]-cart_config_TS[2]),2.));	
	return rtrn_internal_vect;
}

/****************************************************************************
 * Member function to calculate up to quartic contribution in Taylor series 
 * expansion in Diabat ....................................................
****************************************************************************/


double Diabat::TAYLOR_quartic(const std::vector<double>& Q_eqm_Morse_grid_4_any_Diabat,
		PES_catecholate* pesobj) const {
	double diabat_quartic = pesobj->V_eq;
	/* Quadratic contribution */
	for (size_t i=0; i<this->final_Morse_fact_quad.size(); i++)	{
		diabat_quartic += this->final_Morse_fact_quad[i][2]*Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_quad[i][0])-1]*
			Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_quad[i][1])-1];
	}
	/* Cubic contribution */
	for (size_t i=0; i<this->final_Morse_fact_cub.size(); i++)	{
		diabat_quartic += this->final_Morse_fact_cub[i][3]*Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_cub[i][0])-1]*
			Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_cub[i][1])-1]*
			Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_cub[i][2])-1];
	}
	/* Quartic contribution */
	for (size_t i=0; i<this->final_Morse_fact_quartic.size(); i++)	{
		diabat_quartic += this->final_Morse_fact_quartic[i][4]*Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_quartic[i][0])-1]*
			Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_quartic[i][1])-1]*Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_quartic[i][2])-1]*
			Q_eqm_Morse_grid_4_any_Diabat[static_cast<int>(this->final_Morse_fact_quartic[i][3])-1];
	}
	return diabat_quartic;
}


/******************************************************************
 * Member function to calculate beyond quartic corrections in the 
 * Diabat Taylor expansion (Both diagonal & off-diagonal) ........
*******************************************************************/
double Diabat::TAYLOR_octic_tot_crr(const std::vector<double>& Q_eqm_Morse_grid_4_any_Diabat) const {
	double e3_dg = 0.;
	double e33_dg = 0.;
	double e3_e33_od = 0.;
	double morse_quintic_33 = -39998.00303255592;
	double morse_sextic_33 = -213448.10441492498;
	double fit_quintic3 = -0.107358;
	double fit_sextic3 = -0.00975234;
	double fit_septic3 = -0.00862175;
	double fit_octic3 = 0.00209774;
	double a = Q_eqm_Morse_grid_4_any_Diabat[2];
	double b = Q_eqm_Morse_grid_4_any_Diabat[32];
	if (a != 0. && b == 0.)	{
		e3_dg += (fit_quintic3*(pow(a,5)))+(fit_sextic3*(pow(a,6)))+(fit_septic3*(pow(a,7)))+(fit_octic3*(pow(a,8)));
	}
	else if (a == 0. && b != 0.)	{
		e33_dg += (morse_quintic_33/120)*(pow(b,5))+(morse_sextic_33/720)*(pow(b,6));
	}
	else if (a != 0. && b!= 0.)	{
		e3_dg += (fit_quintic3*pow(a,5))+(fit_sextic3*pow(a,6))+(fit_septic3*pow(a,7))+(fit_octic3*pow(a,8));
		e33_dg += (morse_quintic_33/120)*pow(b,5)+(morse_sextic_33/720)*pow(b,6);
	}
	for (size_t i=0; i<this->QS_CRR_OFF_DIAG.size(); i++)	{
//		int l1 = static_cast<int>(this->QS_CRR_OFF_DIAG[i][0]);
//		int l2 = static_cast<int>(this->QS_CRR_OFF_DIAG[i][1]);
		double fact1 = static_cast<double>(std::tgamma(static_cast<int>(this->QS_CRR_OFF_DIAG[i][0])+1));
		double fact2 = static_cast<double>(std::tgamma(static_cast<int>(this->QS_CRR_OFF_DIAG[i][1])+1));
		e3_e33_od += (1/(fact1*fact2))*(this->QS_CRR_OFF_DIAG[i][2])*(pow(a,static_cast<int>(this->QS_CRR_OFF_DIAG[i][0])))*
			(pow(b,static_cast<int>(this->QS_CRR_OFF_DIAG[i][1])));
	}
	return e3_dg+e33_dg+e3_e33_od;
}

/*********************************************************************
 * Member function to calculate total V11 & V22 including asymptotic 
 * correction terms ................................................
**********************************************************************/

std::vector<double> Diabat::Get_tot_V11_and_V22(const std::vector<double>& Q_TS_scaled_cur,
		PES_catecholate* pesobj) const {
	// getting Q_TS mirror grid //
	std::vector<double> Q_TS_mirror_grid_cur = this->Get_Q_TS_mirror_grid(Q_TS_scaled_cur,pesobj);
	std::vector<std::vector<double>> all_internal_vect = this->Get_all_internal_params(Q_TS_scaled_cur,pesobj);
	std::vector<double> Q_EQM_grid = Get_Q_EQM_4_TS_or_Mirror(Q_TS_scaled_cur,pesobj);
	std::vector<double> Q_EQM_MORSE_grid = this->Get_Q_eqm_morse_grid(Q_EQM_grid,pesobj);  
	std::vector<double> Q_EQM_mirror_grid = Get_Q_EQM_4_TS_or_Mirror(Q_TS_mirror_grid_cur,pesobj);
	std::vector<double> Q_EQM_mirror_MORSE_grid = this->Get_Q_eqm_morse_grid(Q_EQM_mirror_grid,pesobj);
	double r_O4H9 = all_internal_vect[0][0];	
	double r_O1H13 = all_internal_vect[1][0];
	double r_O4H10 = all_internal_vect[1][1];
	double aO1C2C3 = all_internal_vect[2][0];
	double aO4C3C2 = all_internal_vect[2][1];
	double r_C2H9 = all_internal_vect[3][0];
	double r_C3H9 = all_internal_vect[3][1];
	double aC2O1H9 = all_internal_vect[3][2];
	double aC3O4H9 = all_internal_vect[3][3];
	double r_O1H9 = all_internal_vect[4][0];
	/* -------------------------------------------*/
	//----------------- V11 ---------------------//
	// Quartic contribution to V11 //
	double V11_quartic = this->TAYLOR_quartic(Q_EQM_MORSE_grid,pesobj); 
	double V11_octic = this->TAYLOR_octic_tot_crr(Q_EQM_MORSE_grid);
	// Non-Bonded O--H (O4-H9) correction V11 //
	double V11_asymp_OH_crr = (60000./2)*(1+tanh((-r_O4H9+0.73)/0.12));
	// O1H13 correction to V11 //
	double V11_asymp_O1H13_crr = 19000./(1+(exp((r_O1H13-2.48)/(2.48/90)))); 
	// O-O correction to V11 //
	double V11_OO_swtch_cur = Q_TS_scaled_cur[0]*this-> mean_slope_OO + this->intercept_OO;
	double V11_asymp_O1O4_crr = 120000./((exp((aO4C3C2-100.8)/1.90))+1);
	double V11_asymp_O1O4_crr_1 = 1/((exp((Q_TS_scaled_cur[9]-V11_OO_swtch_cur)/1.1))+1);
	double V11_asymp_O1O4_crr_2 = V11_asymp_O1O4_crr*V11_asymp_O1O4_crr_1;
	// aC3O4H9 correction to V11 //
	double V11_asymp_aC3O4H9_crr = 40.*((exp((-106.0+aC3O4H9)/5.))+1)*((tanh((-113.+aC3O4H9)/4.))+1);
	/* -------------------------------------------*/
	//----------------- V22 ---------------------//
	// Quartic contribution to V22 //
	double V22_quartic = this->TAYLOR_quartic(Q_EQM_mirror_MORSE_grid,pesobj);
	double V22_octic = this->TAYLOR_octic_tot_crr(Q_EQM_mirror_MORSE_grid);	
	// Non-bonded O--H (O1-H9) correction V22 //
	double V22_asymp_OH_crr = (60000./2)*(1+tanh((-r_O1H9+0.73)/0.12));
	// O1H13 correction to V22 //
	double V11_asymp_O1H13_mirror_crr = 19000./(1+(exp((r_O4H10-2.48)/(2.48/90))));
	// O-O correction to V22 //
	double V22_OO_swtch_cur = Q_TS_mirror_grid_cur[0]*this-> mean_slope_OO + this->intercept_OO;
	double V11_asymp_O1O4_mirror_crr = 120000./((exp((aO1C2C3-100.8)/1.90))+1); 
	double V11_asymp_O1O4_mirror_crr_1 = 1/((exp((Q_TS_scaled_cur[9]-V22_OO_swtch_cur)/1.1))+1); 
	double V11_asymp_O1O4_mirror_crr_2 = V11_asymp_O1O4_mirror_crr*V11_asymp_O1O4_mirror_crr_1;
	// aC2O1H9 correction to V22 //
	double V22_asymp_aC2O1H9_crr = 40.*((exp((-106.0+aC2O1H9)/5.))+1)*((tanh((-113.+aC2O1H9)/4.))+1);
	/* -------------------------------------------*/
	// asymp C2H9/C3H9 correction to V11 and V22 //
	// if Q1 is -ve current grid is used, else corresponding mirror grid is used //
	double asymp_C2H9_swtch_cur;
	if (Q_TS_scaled_cur[0] < 0.)	{
		asymp_C2H9_swtch_cur = Q_TS_scaled_cur[0]*this-> mean_slope_C2H9 + this-> intercept_C2H9; 
	}
	else if (Q_TS_scaled_cur[0] > 0.)	{
		asymp_C2H9_swtch_cur = Q_TS_mirror_grid_cur[0]*this-> mean_slope_C2H9 + this-> intercept_C2H9;
	}
	else if (Q_TS_scaled_cur[0] == 0.)	{
		asymp_C2H9_swtch_cur = 1.128;
	}
	double V11_asymp_C2H9_crr = (6.00e+06)*(exp(1.55-r_C2H9)/0.25)*(1+tanh((asymp_C2H9_swtch_cur-r_C2H9)/0.066));
	double V22_asymp_C3H9_crr = (6.00e+06)*(exp(1.55-r_C3H9)/0.25)*(1+tanh((asymp_C2H9_swtch_cur-r_C3H9)/0.066));
	/* -------------------------------------------*/
	double V11_tot = V11_quartic+V11_octic+V11_asymp_OH_crr+V11_asymp_O1H13_crr
		+V11_asymp_O1O4_crr_2+V11_asymp_C2H9_crr+V11_asymp_aC3O4H9_crr;
	double V22_tot = V22_quartic+V22_octic+V22_asymp_OH_crr+V11_asymp_O1H13_mirror_crr
		+V11_asymp_O1O4_mirror_crr_2+V22_asymp_C3H9_crr+V22_asymp_aC2O1H9_crr;
	std::vector<double> V11_V22 = {V11_tot,V22_tot};
	return V11_V22;
}



/**********************************************************************
 * Function for generating Full potential (Adiabat) corresponding to a 
 * multidimensional grid ...
 * Diabat (V11 & V22) values and V12_sqr are given as arguments
***********************************************************************/

double Get_Adiabat(std::vector<double>& V11_V22, double& V12_sqr)	{
	double E_adiabat;
	double V11_V22_sqr_diff = 0.25*(pow(V11_V22[0] - V11_V22[1],2.));
	if (V11_V22_sqr_diff < 5.0)	{
		if (V12_sqr <0. && std::abs(V12_sqr) < 40.)	{
			V12_sqr = 0.;
			E_adiabat = (0.5*(V11_V22[0]+V11_V22[1]) - sqrt(V11_V22_sqr_diff+V12_sqr));
		}
		else {
			E_adiabat = (0.5*(V11_V22[0]+V11_V22[1]) - sqrt(V11_V22_sqr_diff+V12_sqr));
		}
	}
	else if (V11_V22_sqr_diff > 5.0 && (V11_V22_sqr_diff+V12_sqr)< 0.)	{
		E_adiabat = 0.5*(V11_V22[0]+V11_V22[1]);
	}
	else {
			E_adiabat = (0.5*(V11_V22[0]+V11_V22[1]) - sqrt(V11_V22_sqr_diff+V12_sqr));
	}
	return E_adiabat;
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

