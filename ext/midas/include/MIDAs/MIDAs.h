#include <fstream>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <functional>
#include <vector>
#include <ctime>
#include <math.h>
#include <cstdlib>
#include <sys/time.h> 
#include <unistd.h> 
#include <complex>
//#include <gsl/gsl_linalg.h>
//#include <gsl/gsl_math.h> 
//#include <gsl/gsl_complex.h> 
//#include <gsl/gsl_complex_math.h> 
//#include <gsl/gsl_eigen.h> 
using namespace std;

#ifndef MIDAs_H
#define MIDAs_H

class MIDAs
{
 public:
       
  //Arguments are:(charge state, fine-grained resolution, coarse-grained resolution, min_probability, method used to computed isotopic distribution)
  MIDAs(int charge_state=0, double fine_resolution=0.01, double coarse_resolution=1.0,double cgid_min_prob=1e-150,int flag_method=1);
  
  //Initializing Molecular Formula,
  //Parameters(molecular formula (chemical formula,peptide/protein), modification of cysteine, chemical group n-terminal,chemical group c-terminal,user specified mw and probability)
  void Initialize_Elemental_Composition(string substance, string mc="C", string nt="", string ct="",int flag_substance=2,vector<struct Composition>  user_element_composition = vector<struct Composition> ());

  //Computes FGID
  vector<struct Isotopic_Distribution> Fine_Grained_Isotopic_Distribution();

  //Computes CGID
  vector<struct Isotopic_Distribution> Coarse_Grained_Isotopic_Distribution(); 

  //Return elements isotopic
  vector< vector<struct Composition> > ELEMENTAL_ISOTOPIC_COMPOSITION();

  //Computing elemental composition from computed isotopic distribution. Under Investigation
  // vector<struct Composition> ELEMENTAL_COMPOSITION_ID(vector<struct Isotopic_Distribution> &,vector<int> &,int,double);
  //void  MATRIX_INV_ID_LS_AVE(vector<struct Isotopic_Distribution> &,vector<int> &,int);
  //void  MATRIX_INV_ID_LS_RATIO(vector<struct Isotopic_Distribution> &,vector<int> &,int);
  
  //Printing FGID to specified file
  void Print_Fine_Structure_Isotopic_Distribution(char *);
  //Printing CGID to specified file
  void Print_Coarse_Structure_Isotopic_Distribution(char *);

  void Elements_Variance_Mass(vector<double> &);
  void Elements_Average_Mass(vector<double> &);
  void Elements_Raw_Moments(vector< vector<double> > &,int);
  void Elements_Cumulants(vector< vector<double> > &,int);
  double monoisotopic_mass();        //returns monoisotopic mass
  double average_molecular_mass();   //returns theoretical average mass
  double variance_molecular_mass();   //returns theoretical average mass
  double average_molecular_mass_distribution(vector<struct Isotopic_Distribution> &);  //returns distribution average mass 
  double variance_molecular_mass_distribution(vector<struct Isotopic_Distribution> &); //returns distribution variance
  void compute_distribution_cumulant_moments(vector<struct Isotopic_Distribution> &,vector<double> &,int);
  void compute_distribution_central_moments(vector<struct Isotopic_Distribution> &,vector<double> &,int);
  double fgid_time();    //returns time of FGID
  double cgid_time();    //returns time of CGID


 private:

  //Computing elemental composition. Used to check computed elemental composition
  vector<int> check_elemental_composition(vector< vector<double> > &,vector<double> &,vector<double> &,vector<int> &,vector<int> &,double&,double,int); 

  //Returns monoisotopic mass given a elemental composition
  double monoisotopic_molecular_mass(vector< vector<struct Composition> > &); 

  /*Elemental Composition*/
  void add_charge_state(vector< vector<struct Composition> > &,int);
  void sort_element_composition(vector< vector<struct Composition> > &);

  /*Adjusting resolution*/
  void set_resolution(double);

  /*Polynomial Expansion for Coarse-Grained Isotopic Distriubtion*/
  vector<struct Polynomial> multiply_polynomial(vector< vector<struct Composition> > &);
  void multiplication_operation( vector<struct Composition> &, vector<struct Polynomial> & );
  void multiplication_final_operation( vector<struct Polynomial> &, vector<struct Polynomial> & );
  vector<struct Polynomial> Merge_Coarse_Polynomial(vector<struct Polynomial>);
 
  /*Polynomial Expansion for Fine-Grained Isotopic Distribution*/
  vector<struct Polynomial> multiply_fine_polynomial(vector< vector<struct Composition> > &);
  void multiplication_fine_final_operation( vector<struct Polynomial> &, vector<struct Polynomial> & );
  vector<struct Polynomial> Merge_Fine_Polynomial(vector<struct Polynomial>);
 
  /*Fourier Transform Functions*/ 
  void FT_Coarse_Grained_ID(vector< vector<struct Composition> > &,vector<struct Isotopic_Distribution> &,double); 
  void FT_Fine_Grained_ID(vector< vector<struct Composition> > &,vector<struct Polynomial> &,double); 
  void FT_JFC_ID(vector< vector<struct Composition> > &);
  void four1(double *,int, int);
  double ft_average_molecular_mass(vector< vector<struct Composition> > &,double resolution=0.0);   //returns theoretical average mass for FT
  double ft_variance_molecular_mass(vector< vector<struct Composition> > &,double resolution=0.0);   //returns theoretical average mass for FT
 
  //Computes Factorial
  double FACTR_LN(int);
  double GAMMA_LN(double xx);

  //Global Variables
  unsigned int VECTOR_CGID_SIZE,VECTOR_FGID_SIZE;
  int CHARGE_STATE,FINE_GRID,FLAG_METHOD;
  double MIN_PROB;
  double MW_RESOLUTION;
  double FINE_RESOLUTION;
  double FINE_MIN_PROB;
  double MASS_HYDROGEN_ION;
  double MASS_ELECTRON;
  double COARSE_RESOLUTION;
  double MERGE_COARSE_RESOLUTION;
  double MERGE_FINE_RESOLUTION;
  double FGID_TIME,CGID_TIME;
  double MONOISOTOPIC_MASS;  

  double *MASS_FREQUENCY_DOMAIN,*MASS_FREQUENCY_DOMAIN_DP; 
  double **ELEMENT_FUNCTION;

  struct Complex *FFT_DATA;

  vector<struct Polynomial> CGID_Polynomial,FGID_Polynomial;
  vector< vector<struct Composition> > Element_Composition;
  vector<struct Composition>  User_Element_Composition;
  vector<struct Isotopic_Distribution> FINE_GRAINED_ISOTOPIC_DISTRIBUTION;
  vector<struct Isotopic_Distribution> COARSE_GRAINED_ISOTOPIC_DISTRIBUTION;
};

//Store the molecule elemental composition
struct Composition
{
  char element[3];
  double power;
  double prob;
  double log_prob;
  double mw;
  double nucleon;
  double ave_mw;
  int atoms;
};
struct by_composiont_mw
{
   inline bool operator ()(const struct Composition &e,const struct Composition &f)
  const {return  f.mw < e.mw;} 
};

//Store the multiplied polynomial
struct Polynomial
{
  double power;
  double prob;
  inline bool operator < (const struct Polynomial &e) const { return power < e.power;}  //increasing order of molecular weight 
};
struct by_prob 
{ 
  inline bool operator ()(const struct Polynomial &e,const struct Polynomial &f)
  const {return  f.prob < e.prob;}
};

//Store final isotopic distribution
struct Isotopic_Distribution
{
  double mw;
  double prob;
  inline bool operator < (const struct Isotopic_Distribution &e) const { return mw < e.mw;}  //increasing order of molecular weight 
};
 
//Used for FFT
struct Complex
{
  double real;
  double imaginary;
};



#endif
