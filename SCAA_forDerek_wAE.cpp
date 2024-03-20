//Stock Catch at age assessment
//created by Nicholas Fisch to be used in simulation study, CKMR data to be added later

//Dec 2022

#include <TMB.hpp>

//this namespace is used for printing during iterations using CppAD::PrintFor
 namespace CppAD {
  void PrintFor(const char* before, const double& var) {  }
}

template<class Type>
Type objective_function<Type>::operator() ()
{

  //TABLE OF CONTENTS
  //I. DATA INPUTS
  //II. PARAMETER DECLARATION
  //III. SETTING INTITAL VALUES FOR DATA AND BACK TRASFORMING PARAMETERS ("Pre-function" section)
  //IV. FUNCTIONS AND CALCULATIONS
  //    1. Selectivity
  //    2. Mortality
  //    3. Population
  //    4. Expected Catch and Catch at Age and index
  //    5. Objective Functions
  //V. REPORT SECTION

  //I. DATA INPUTS
  //---------------------------------------------------------------------------------------------
  DATA_INTEGER(fyear); //First year
  DATA_INTEGER(lyear); //Last year
  DATA_INTEGER(fage); //First age
  DATA_INTEGER(lage); //Last age

  DATA_VECTOR(years); //Vector of years
  DATA_VECTOR(ages); //Vector of ages

  DATA_VECTOR(obs_harv);              //Observed harvest
  DATA_VECTOR(obs_index);             //Observed index
  DATA_MATRIX(obs_fishery_comp);      //Observed fishery Composition
  DATA_VECTOR(SS_fishery);            //Sample size for fishery compositions

  DATA_VECTOR(Mat);                 //Maturity at age vector 
  DATA_VECTOR(Laa);                 //Length at age vector
  DATA_VECTOR(Waa);                 //Weight at age vector

  DATA_INTEGER(Lamda_Harvest);      //Whether to use data sources or not
  DATA_INTEGER(Lamda_Comp);    
  DATA_INTEGER(Lamda_Index);    

  DATA_MATRIX(AE_mat);              //Ageaing Error Matrix

  //II. PARAMETER DECLARATION
  //---------------------------------------------------------------------------------------------
  PARAMETER(log_M);          //M scalar for Lorenzen M
  PARAMETER(log_q);          //Fishery Catchability for index
  
  PARAMETER_VECTOR(log_recruit_devs_init); //Log scale recruitment deviations for initial abundance
  PARAMETER_VECTOR(log_recruit_devs);      //Log scale recruitment deviations
  PARAMETER(steepness);                    //Steepness for recruitment
  PARAMETER(log_R0);                       //Unfished recruitment
  PARAMETER(log_sigma_rec);                //Log sd for recruitment

  PARAMETER(log_sd_catch);           //Log sd for catch
  PARAMETER(log_sd_index);           //Log sd for index

//  PARAMETER(Sel_logis_k);         //logistic selectivity slope parameter 
//  PARAMETER(Sel_logis_midpt);     //logistic selectivity midpoint parameter
  PARAMETER(B1);                     //Double normal sel pars
  PARAMETER(B2);     
  PARAMETER(B3);         
  PARAMETER(B4);     
  PARAMETER(B5);         
  PARAMETER(B6);     

  PARAMETER_VECTOR(log_fint);       //Log scale fishing intensities (fully selected fishing mortalities)
  //---------------------------------------------------------------------------------------------

//DERIVED QUANTITIES
 
//back-transform log-scale params
  Type NatMort=exp(log_M);
  Type q=exp(log_q);
  Type R0=exp(log_R0);
  Type sd_rec=exp(log_sigma_rec);                //sigma for recruit deviations
  Type sd_catch=exp(log_sd_catch);           //sigma for fishery catch
  Type sd_index=exp(log_sd_index);               //cv for fishery CPUE

  vector<Type> fishery_sel(ages.size());                   //Fishery selectivity (fage,lage)

  vector<Type> Maa(ages.size());                      //Natural mortality vector
  matrix<Type> M(years.size(),ages.size());                //Instantaneous natural mortality matrix
  matrix<Type> Sel(years.size(),ages.size());              //Selectivity
  matrix<Type> F(years.size(),ages.size());                //Instantaneous fishing mortality matrix
  matrix<Type> Z(years.size(),ages.size());                //Total mortality matrix
  matrix<Type> S(years.size(),ages.size());                //Annual survival matrix
  matrix<Type> A(years.size(),ages.size());                //Annual mortality matrix

  matrix<Type> N(years.size()+1,ages.size());	           //Predicted abundance at age, +1 is to incorporate final catch 
  matrix<Type> Nbar(years.size(),ages.size());	           //Nbar

  matrix<Type> pred_caa(years.size(),ages.size());              //Predicted fishery catch at age
  matrix<Type> pred_fishery_comp(years.size(),ages.size());     //Predicted true age composition from the fishery catch
  matrix<Type> pred_fishery_comp_wAE(years.size(),ages.size()); //Predicted age composition from the fishery catch, with ageing error applied
  vector<Type> pred_harv(years.size());                         //Predicted Harvest weight (kg), (fyear,lyear) 

  vector<Type> pred_index(years.size());                    //Predicted index, (fyear,lyear)
  
  vector<Type> log_rec_devs(years.size()+ages.size());    //rec devs vector that adds in the current year, (fage,lyear+lage)

  vector<Type> spbiomass(years.size()+1);                  //Spawning biomass, (fyear,lyear+1)
  vector<Type> N0_age(ages.size());                        //Unfished numbers at, (fage,lage)
  vector<Type> lxo(ages.size());                           //Unfished numbers at age, (fage,lage)
  Type SSB0;                                               //Unfished spawning biomass
       
//Double normal stuff   
  Type peak2; 
  Type t1; 
  Type t2; 
  vector<Type> j1(ages.size());
  vector<Type> j2(ages.size());
  vector<Type> asc(ages.size());
  vector<Type> dsc(ages.size());
  
  //Likelihood components
  Type L1;
  Type L2;
  Type L3;
  Type L4;
 
  int i; 
  int j; 
  
  Type NLL = 0;
  Type NPRAND = 0;
  Type JNLL = 0;
  
//////////////////////////////////////
///////////////////////////
//Actual Model
///////////////////////////
/////////////////////////////////////

///////////////////////////
//FISHERY SELECTIVITY
///////////////////////////

//Logistic selectivity
//  for(j=0;j<=lage;j++){
//   fishery_sel(j)=1/(1+exp(-log(19)*(Laa(j)-exp(Sel_logis_midpt))/exp(Sel_logis_k)));   //Mean length at age Logistic Selectivity, stock synthesis version
//   fishery_sel(j)=1/(1+exp(-log(19)*(int(j)-exp(Sel_logis_midpt))/exp(Sel_logis_k)));   //Age-based Logistic Selectivity, stock synthesis version
//  }
  
  //Double Normal
  Type Amin=fage;
  Type Amax=lage;

  for(j=0;j<=ages.size()-1;j++){
   peak2=B1+1+((0.99*Amax-B1-1)/(1+exp(-B2)));
   t1=exp(-pow(Amin-B1,2)/exp(B3));
   t2=exp(-pow(Amax-peak2,2)/exp(B4));
  
   j1(j)=pow((1+exp(-20*((int(j)-B1)/(1+fabs(int(j)-B1))))),-1);
   j2(j)=pow((1+exp(-20*((int(j)-peak2)/(1+fabs(int(j)-peak2))))),-1);

   //Long form
   asc(j)=pow(1+exp(-B5),-1)+(1-pow(1+exp(-B5),-1))*((exp(-pow(int(j)-B1,2)/exp(B3))-t1)/(1-t1));
   dsc(j)=1+(pow(1+exp(-B6),-1)-1)*((exp(-(int(j)-peak2)/exp(B4))-1)/(t2-1));

   //Actual Sel
   fishery_sel(j)=asc(j)*(1-j1(j))+j1(j)*((1-j2(j))+j2(j)*dsc(j));
  }
  
//////////////////////
//MORTALITY
//////////////////////
 //Lorenzen M or Constant M
//  for(j=0;j<=lage;j++){
//   Maa(j)=NatMort*pow((Lt(j)/(Linf*0.75)),-1);
// }
   Maa.fill(NatMort);

  for(i=0;i<=years.size()-1;i++){
   for(j=0;j<=lage;j++){
    F(i,j) = fishery_sel(j)*exp(log_fint(i)); //setting year and age specific fishing mortality as the product of selectivity and year specific fishing intensity
	Z(i,j) = Maa(j) + F(i,j);                 //Total instantaneous mortality
    S(i,j) = exp(-1.0*Z(i,j));                //Annual survival
    A(i,j)=1.0-S(i,j);                        //Annual mortality
   }
  }
 
//////////////////////
//POPULATION
//////////////////////
  //Unfished Spawning biomass calcs
  lxo(0)=1.0;
  for(j=fage+1;j<=lage;j++){
   lxo(j)=lxo(j-1)*exp(-1*Maa(j-1));   //cumulative survival
  }
  lxo(lage)=lxo(lage)/(1-exp(-1*Maa(lage)));
  N0_age=R0*lxo;
  SSB0=(N0_age*Mat*Waa).sum();   //Numbers at age * Fecundity

  //Filling in log_rec_devs
  for(i=0;i<=years.size()-1;i++){
   log_rec_devs(i)=log_recruit_devs(i);
  }
  log_rec_devs(lyear)=Type(0);
  
  //Abundance at age in the first year
  for(j=fage+1;j<=lage;j++){
//   N(fyear-1,j)=R0*lxo(j);          //set abundance in the fist year f
   N(fyear-1,j)=R0*exp(log_recruit_devs_init(fabs(j-(lage)))-0.5*pow(sd_rec,2))*lxo(j);          //set abundance in the fist year f
  }
//  N(fyear-1,fage)=R0*exp(log_rec_devs(0));     //Filling in initial year recruitment
  N(fyear-1,fage)=R0*exp(log_rec_devs(0)-0.5*pow(sd_rec,2));     //Filling in initial year recruitment
  
  //Spawning Biomass in the first year
  spbiomass(fyear-1)=(vector<Type>(N.row(fyear-1))*Mat*Waa).sum();
  
//Population loop
  for(i=fyear;i<=years.size();i++){  //TMB starts at zero so this is 81 years (or 41), starting at i+1
    for(j=fage+1;j<=lage;j++){
      N(i,j)=N(i-1,j-1)*S(i-1,j-1);
    }
   N(i,lage)+=N(i-1,lage)*S(i-1,lage);   //Plus group
   
  spbiomass(i)=(vector<Type>(N.row(i))*Mat*Waa).sum();

//   N(i,fage)=((4.*steepness*R0*spbiomass(i))/(SSB0*(1.-steepness)+spbiomass(i)*(5.*steepness-1.)))*exp(log_rec_devs(i));  //Recruitment if at age 0
   N(i,fage)=((4.*steepness*R0*spbiomass(i))/(SSB0*(1.-steepness)+spbiomass(i)*(5.*steepness-1.)))*exp(log_rec_devs(i)-0.5*pow(sd_rec,2));  //Recruitment, if at age 1
  }
 
////////////////
//CATCH
////////////////

  pred_caa=F.array()/Z.array()*A.array()*N.block(0,fage,years.size(),lage+1).array();
  Nbar=N.block(0,fage,years.size(),lage+1).array()*((1.0-exp(-1.0*Z.array())))*(1.0/Z.array());
  
  for(i=0;i<=years.size()-1;i++){  //TMB starts at zero   
   pred_harv(i)=(vector<Type>(pred_caa.row(i))*Waa).sum();
   pred_index(i)=q*(vector<Type>(Nbar.row(i))*Waa*fishery_sel).sum();
   pred_fishery_comp.row(i)=pred_caa.row(i)/(pred_caa.row(i).sum());  //calculating predicted catch age composition
  }

  pred_fishery_comp_wAE=pred_fishery_comp*AE_mat;
  
/////////////////////////
//Objective function
/////////////////////////
  Type pi=3.141593;

//LIKELIHOODS
  L1=Type(0);
  //Catch Likelihood
  // Lognormal by SD
  for(i=0;i<=lyear-1;i++){  
   L1 += log(obs_harv(i))+0.5*log(2*pi)+log_sd_catch+pow(log(obs_harv(i))-log(pred_harv(i)),2)/(2*pow(sd_catch,2));
  } 

//Multinomial
  L2=Type(0);
  for(i=0;i<=lyear-1;i++){
   if(SS_fishery(i)>0){   //Only run calcs if there is comp data
    for(j=fage;j<=lage;j++){
     L2 += -1*(SS_fishery(i)*(obs_fishery_comp(i,j)*log(pred_fishery_comp_wAE(i,j))));     //Likelihood for age composition of fishery catch
    }
   }
  }

//Fishery Index
//Lognormal by SD
  L3=Type(0);
  for(i=0;i<=lyear-1;i++){
   L3 += log(obs_index(i))+0.5*log(2*pi)+log_sd_index+pow(log(obs_index(i))-log(pred_index(i)),2)/(2*pow(sd_index,2));
  }

  //Recruitment deviations
   NPRAND = -1*(sum(dnorm(log_recruit_devs,0,sd_rec,true)) + sum(dnorm(log_recruit_devs_init,0,sd_rec,true))); //Recruitment deviations

  NLL = (Lamda_Harvest*L1) + (Lamda_Comp*L2) + (Lamda_Index*L3);
  JNLL=NPRAND+NLL;

/////////////////////////
//Report
////////////////////////////

  REPORT(Maa);
  REPORT(log_M);
  REPORT(q);
  REPORT(R0);
  REPORT(sd_rec);
  REPORT(sd_catch);
  REPORT(sd_index);
  REPORT(log_recruit_devs);
  REPORT(B1);                     //Double normal sel pars
  REPORT(B2);     
  REPORT(B3);         
  REPORT(B4);     
  REPORT(B5);         
  REPORT(B6);     

  REPORT(fishery_sel);
  REPORT(N0_age);
  REPORT(Laa);
  REPORT(lxo);
  REPORT(F);
  REPORT(M);
  REPORT(N);
  REPORT(spbiomass);
  ADREPORT(spbiomass);
  REPORT(obs_harv);
  REPORT(pred_harv);
  REPORT(pred_caa);
  REPORT(pred_index);
  REPORT(obs_index);
  REPORT(pred_fishery_comp);
  REPORT(pred_fishery_comp_wAE);
  REPORT(obs_fishery_comp);

  REPORT(L1);
  REPORT(L2);
  REPORT(L3);
  REPORT(L4);
  REPORT(NLL);
  REPORT(NPRAND);

  return JNLL;

}

