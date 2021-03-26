# if ! defined(_PES_H)
# define _PES_H

# include <iosfwd>
# include <vector>
# include <string>


/***************************************************************
 * class PES_catecholate contains private data members specific
 * to the PES (most of which are constants provided in relevant 
 * header for constants
 * member function V12_sqr calculates coupling value for a given
 * multidimensional grid 
 * *************************************************************/


class PES_catecholate	{
	private:
		double prefact_4_scaled_normal_coords;	// prefactor for generating scale factors for transforming to scaled normal coordinates 
		double V_eq;							// EQM energy in Hartree
		double morse_a;
		std::vector<int> Mode_flag_TS;			// A1 modes are flagged as 1 and B2,B1,A2 modes are flagged as -1
		std::vector<double> REF_TS;				// COM removed TS cartesians 
		std::vector<double> REF_EQM;			// COM removed EQM cartesians
		std::vector<double> l_TS;				// TS l-matrix as obtained from Gaussian 09
		std::vector<double> EQM_rot_mat;
		std::vector<double> l_eqm;				// EQM l-matrix (mass-unweighted & renormalized)
		std::vector<double> mass_sqrt;
		std::vector<double> ratio_TS;
		std::vector<std::vector<double>> Q_K_TS_scaled;
		std::vector<double> freq_EQM; 
		std::vector<std::vector<std::vector<double>>> B_fit_vect;
		std::vector<double> alpha_vect;
	public:
		// constructor for PES_catecholate 
		PES_catecholate(double prefact_4_scaled_normal_coords,double V_eq,double morse_a,
				std::vector<int>& Mode_flag_TS,std::vector<double>& REF_TS,
				std::vector<double>& REF_EQM,std::vector<double>& l_TS,
				std::vector<double>& EQM_rot_mat,std::vector<double>& l_eqm,
				std::vector<double>& mass_sqrt,std::vector<double>& ratio_TS,
				std::vector<std::vector<double>>& Q_K_TS_scaled,std::vector<double>& freq_EQM,
				std::vector<std::vector<std::vector<double>>>& B_fit_vect,
				const std::vector<double>& alpha_vect);
		std::vector<double> TS_scaled_2_Gauss(const std::vector<double>& Q_TS_scl) const;
		std::vector<double> Eqm_rnrm_2_scaled(const std::vector<double>& Q_eq_rnrm) const;
		double V12_sqr(const std::vector<double>& Q_TS_scaled_cur) const;
		/* Diabat class friend of PES_catecholate */
		friend class Diabat;
};


/***************************************************************
 * Class Diabat contains numerical derivatives for generating 
 * Taylor expansion .. and some constants for asymptotic 
 * correction ..
 * Member function Get_tot_V11_and_V22 calculates total V11 and 
 * V22 contribution for a given multidimensional grid 
***************************************************************/

class Diabat	{
	private:
		std::vector<std::vector<double>> final_Morse_fact_quad;
		std::vector<std::vector<double>> final_Morse_fact_cub;
		std::vector<std::vector<double>> final_Morse_fact_quartic;
		std::vector<std::vector<double>> QS_CRR_OFF_DIAG;
		double mean_slope_OO;
		double intercept_OO;
		double mean_slope_C2H9;
		double intercept_C2H9;
		public:
			Diabat(std::vector<std::vector<double>>& final_Morse_fact_quad,std::vector<std::vector<double>>& final_Morse_fact_cub,
					std::vector<std::vector<double>>& final_Morse_fact_quartic,std::vector<std::vector<double>>& QS_CRR_OFF_DIAG,
					double mean_slope_OO,double intercept_OO,double mean_slope_C2H9,double intercept_C2H9);
			std::vector<double> Get_Q_TS_mirror_grid(const std::vector<double>& Q_TS_scaled_cur,PES_catecholate* pesobj) const;
			std::vector<double> Get_Q_EQM_4_TS_or_Mirror(const std::vector<double>& Q_TS_scaled_cur,PES_catecholate* pesobj) const;
			std::vector<double> Get_Q_eqm_morse_grid(const std::vector<double>& Q_eqm_grid,PES_catecholate* pesobj) const;
			std::vector<std::vector<double>> Get_all_internal_params(const std::vector<double>& Q_TS_scaled_cur,PES_catecholate* pesobj) const;
			double TAYLOR_quartic(const std::vector<double>& Q_eqm_Morse_grid_4_any_Diabat,PES_catecholate* pesobj) const;	
			double TAYLOR_octic_tot_crr(const std::vector<double>& Q_eqm_Morse_grid_4_any_Diabat) const;
			std::vector<double> Get_tot_V11_and_V22(const std::vector<double>& Q_TS_scaled_cur,PES_catecholate* pesobj) const;
};



/***************************************************************
 * Function for generating Full potential (Adiabat) corresponding
 * to a multidimensional grid ..
 * Diabat (V11 & V22) values and V12_sqr are given as arguments
***************************************************************/

double Get_Adiabat(std::vector<double>& V11_V22,double& V12_sqr);

/***************************************************
 * Function for generating line number vector from 
 * number of lines ..
 * *************************************************/

std::vector<int> Get_lineno_vect(int NUMLINES);



# endif
