//FILE CREATED ON 08/13/09 BY TRAVIS BRENDEN FOR ADMB 2009 SHORT COURSE
//MODIFIED BY KYLE MOLTON & TRAVIS BRENDEN FOR WHITEFISH MSE SIMULATION PROJECT
//LAST MODIFIED 1/15/2013 By YANG
//SIMPLE SCAA MODEL
GLOBALS_SECTION
  #include <time.h>
  #include <admodel.h>
  //int seed=45456; //fixed random seed
  int seed=(unsigned)time(0); //using system time as random seed
  //the pesudo code I used from
  //http://commons.apache.org/math/apidocs/org/apache/commons/math/stat/descriptive/rank/Percentile.html
  //this match to type=6 in quantile function in R, which is also used by SPSS and minitab,
  // but R by defult use type=7, and R also implemented total 9 types
  //require the index of data start from 1

  double quantile(dvector data, double p){
    int n=data.size();
    dvector sortdata(1,n); 
    sortdata=sort(data);
    //cout<<sortdata<<endl<<endl;

    if(n==1) return sortdata(1);
    double pos=p*(n+1); // estimated percentile position
    double d=pos - floor(pos); //difference
    //cout<<pos<<" "<< d<<endl;

    if(pos<1) return sortdata(1);
    else if(pos>=n) return sortdata(n);
    else{
      double lower=sortdata(floor(pos));
      double upper=sortdata(floor(pos)+1);
      return (lower+ d*(upper-lower));
    }
  }

  double quantile(dvar_vector data, double p){
      return quantile(value( data),p);
    
  }
DATA_SECTION
  int fyear
  !!fyear=1;
  init_int Pops
  init_matrix Mix(1,Pops,1,Pops)
  init_int lyear                                    //Number of years
  init_int fage
  init_int lage                                     //Number of ages
  init_number nat_mort
  init_number tot_targ_Z
  init_vector wt_age(fage,lage)                     //weight at age
  init_vector sel_age(fage,lage)                    //not being used, comment by liu
  init_vector mat_age(fage,lage)        //maturity at age
  init_number propFemale                 //female proportion
  init_number IntP1
  init_number IntP2
  init_number IntP3
  init_number IntP4 
  init_number RecP1
  init_number RecP2
  init_number RecP3
  init_number RecP4       
  init_matrix rc_catch1(fyear,lyear,fage,lage)      //recreational catch
  init_matrix rc_catch2(fyear,lyear,fage,lage)      //recreational catch
  init_matrix rc_catch3(fyear,lyear,fage,lage)      //recreational catch
  init_matrix rc_catch4(fyear,lyear,fage,lage)      //recreational catch
  init_vector rc_eff1(fyear,lyear)		              //observed recreational effort
  init_vector rc_eff2(fyear,lyear)
  init_vector rc_eff3(fyear,lyear)
  init_vector rc_eff4(fyear,lyear)
  init_int lastNumYrMort  //how many last yrs mortality as mean mortality for TAC
  init_int lastNumYrRec    //the last 5yrs first age as mean recruitment
  init_number stdObsCatch  //lognormal std error for observed catch
  init_int effSampleSize  //effective sample size
  init_vector test(1,3)				//test vector for ensuring data are read in correctly
  
  int i                                         //index for year loop
  int j                                         //index for age loop
  matrix M(fyear,lyear,fage,lage)               //Natural mortality matrix
  vector rc_cch1(fyear,lyear)		        //observed recreational catch
  vector rc_cch2(fyear,lyear)
  vector rc_cch3(fyear,lyear)
  vector rc_cch4(fyear,lyear)
  matrix rc_cps1(fyear,lyear,fage,lage)          //age composition from the recreational catch
  matrix rc_cps2(fyear,lyear,fage,lage)
  matrix rc_cps3(fyear,lyear,fage,lage)
  matrix rc_cps4(fyear,lyear,fage,lage)
  matrix obsResampleCatch1(fyear,lyear,fage,lage)  //observed catch based on input actual age composition and effective sample size
  matrix obsResampleCatch2(fyear,lyear,fage,lage)
  matrix obsResampleCatch3(fyear,lyear,fage,lage)
  matrix obsResampleCatch4(fyear,lyear,fage,lage)
  vector new_rc_cch1(fyear,lyear)   //observed annual total catch after add in obs.error
  vector new_rc_cch2(fyear,lyear)
  vector new_rc_cch3(fyear,lyear)
  vector new_rc_cch4(fyear,lyear)
  matrix new_rc_catch1(fyear,lyear,fage,lage)  //obs catch at age based on observed total catch times new age composition
  matrix new_rc_catch2(fyear,lyear,fage,lage)
  matrix new_rc_catch3(fyear,lyear,fage,lage)                
  matrix new_rc_catch4(fyear,lyear,fage,lage)
  matrix new_rc_cps1(fyear,lyear,fage,lage)    //age composition from resampling readin catch
  matrix new_rc_cps2(fyear,lyear,fage,lage)
  matrix new_rc_cps3(fyear,lyear,fage,lage)
  matrix new_rc_cps4(fyear,lyear,fage,lage)
  
  
  number mre1
  number mre2
  number mre3
  number mre4
  number mare1
  number mare2
  number mare3
  number mare4

  //!!cout << test << endl;
  //!!exit(88);


INITIALIZATION_SECTION

  log_qrc -13.41005              //Recreational catchability
  log_rc_s_p1 2.570625525          //Recreational selectivity parameter 1, estimated on a log scale
  rc_s_p2 1.26                   //Recreational selectivity parameter 2
  log_initial1 IntP1
  log_initial2 IntP2
  log_initial3 IntP3
  log_initial4 IntP4
  log_mean_recruit1 RecP1
  log_mean_recruit2 RecP2
  log_mean_recruit3 RecP3
  log_mean_recruit4 RecP4
  log_sigma_rc -1.0




PARAMETER_SECTION

  //Declaration of parameters to be estimated
  init_bounded_number log_qrc(-15.,-10.,1)                 //Recreational catchability coefficient
  init_bounded_number log_rc_s_p1(-3.,5.,1)                //Recreational selectivity parameter 1
  init_bounded_number rc_s_p2(-3.,3.,1)                    //Recreational selectivity parameter 2
  init_bounded_number log_initial1(8.5,19.,1)                 //Numbers of age-2 and older fish in the first year
  init_bounded_number log_initial2(9.,19.,1) 
  init_bounded_number log_initial3(9.,19.,1) 
  init_bounded_number log_initial4(9.,19.,1) 
  init_bounded_dev_vector initial_devs1(fage+1,lage,-5.,5.,1)
  init_bounded_dev_vector initial_devs2(fage+1,lage,-5.,5.,1)
  init_bounded_dev_vector initial_devs3(fage+1,lage,-5.,5.,1)
  init_bounded_dev_vector initial_devs4(fage+1,lage,-5.,5.,1)
  init_bounded_number log_mean_recruit1(9.0,19.,1)           //Average recruitment over the years - log scale
  init_bounded_number log_mean_recruit2(9.0,19.,1) 
  init_bounded_number log_mean_recruit3(9.0,19.,1) 
  init_bounded_number log_mean_recruit4(9.0,19.,1) 
  init_bounded_dev_vector recruit_devs1(fyear,lyear,-3.,3.,1)  //Annual recruitment deviations
  init_bounded_dev_vector recruit_devs2(fyear,lyear,-3.,3.,1)
  init_bounded_dev_vector recruit_devs3(fyear,lyear,-3.,3.,1)
  init_bounded_dev_vector recruit_devs4(fyear,lyear,-3.,3.,1)
  init_bounded_dev_vector log_effort_devs1(fyear,lyear,-20.,20.,1)
  init_bounded_dev_vector log_effort_devs2(fyear,lyear,-20.,20.,1)
  init_bounded_dev_vector log_effort_devs3(fyear,lyear,-20.,20.,1)
  init_bounded_dev_vector log_effort_devs4(fyear,lyear,-20.,20.,1)
  init_bounded_number log_sigma_rc(-5.,5.,1)                 //log sd for rec catch

  
  vector initial1(fage+1,lage)                  //Exponent of age-2 and older fish in the first year
  vector initial2(fage+1,lage)
  vector initial3(fage+1,lage)
  vector initial4(fage+1,lage)
  number qrc                                   //Exponent of recreational catchability
  number rc_s_p1                               //Exponent of recreational selectivity parameter
  vector rc_s_age(fage,lage)                   //Recreational selectivities by age
  matrix rc_s(fyear,lyear,fage,lage)           //Matrix of recreational selectivities
  vector re1(fyear,lyear)
  vector re2(fyear,lyear)
  vector re3(fyear,lyear)
  vector re4(fyear,lyear)
  vector abs_re1(fyear,lyear)
  vector abs_re2(fyear,lyear)
  vector abs_re3(fyear,lyear)
  vector abs_re4(fyear,lyear)
  

  matrix F1(fyear,lyear,fage,lage)              //Instantaneous fishing mortality matrix
  matrix F2(fyear,lyear,fage,lage)
  matrix F3(fyear,lyear,fage,lage)
  matrix F4(fyear,lyear,fage,lage)
  matrix Z1(fyear,lyear,fage,lage)              //Total mortality matrix
  matrix Z2(fyear,lyear,fage,lage)
  matrix Z3(fyear,lyear,fage,lage)
  matrix Z4(fyear,lyear,fage,lage)
  matrix S1(fyear,lyear,fage,lage)              //Annual survival matrix
  matrix S2(fyear,lyear,fage,lage)
  matrix S3(fyear,lyear,fage,lage)
  matrix S4(fyear,lyear,fage,lage)
  matrix n1(fyear,lyear,fage,lage)	       //Predicted abundance at age
  matrix n2(fyear,lyear,fage,lage)
  matrix n3(fyear,lyear,fage,lage)
  matrix n4(fyear,lyear,fage,lage)
  matrix ssb1(fyear,lyear,fage,lage)	       //Predicted abundance weight(~ssb) at age
  matrix ssb2(fyear,lyear,fage,lage)
  matrix ssb3(fyear,lyear,fage,lage)
  matrix ssb4(fyear,lyear,fage,lage)
  vector Total_ssb1(fyear,lyear)                        //Predicted total abundance weight(~ssb)
  vector Total_ssb2(fyear,lyear)
  vector Total_ssb3(fyear,lyear)
  vector Total_ssb4(fyear,lyear)
  vector N1(fyear,lyear)                        //Predicted total abundance
  vector N2(fyear,lyear)
  vector N3(fyear,lyear)
  vector N4(fyear,lyear)
  matrix n_rc1(fyear,lyear,fage,lage)	       //Predicted recreational catch at age
  matrix n_rc2(fyear,lyear,fage,lage)
  matrix n_rc3(fyear,lyear,fage,lage)
  matrix n_rc4(fyear,lyear,fage,lage)
  vector N_rc1(fyear,lyear)                     //Predicted total recreational catch
  vector N_rc2(fyear,lyear)
  vector N_rc3(fyear,lyear)
  vector N_rc4(fyear,lyear)
  matrix n_rc_cps1(fyear,lyear,fage,lage)       //Predicted recreational catch age composition
  matrix n_rc_cps2(fyear,lyear,fage,lage)
  matrix n_rc_cps3(fyear,lyear,fage,lage)
  matrix n_rc_cps4(fyear,lyear,fage,lage)
  //matrix d_M1(fyear,lyear,fage,lage)            //Predicted natural deaths at age
  //matrix d_M2(fyear,lyear,fage,lage)
  //matrix d_M3(fyear,lyear,fage,lage)
  //matrix d_M4(fyear,lyear,fage,lage)
  vector effort_devs1(fyear,lyear)
  vector effort_devs2(fyear,lyear)
  vector effort_devs3(fyear,lyear)
  vector effort_devs4(fyear,lyear)
  vector rc_diff1(fyear,lyear)                  //Differences between observed and predicted recreational catches
  vector rc_diff2(fyear,lyear)
  vector rc_diff3(fyear,lyear)
  vector rc_diff4(fyear,lyear)
  number log_sigma_effort
  number variance_ratio


  number L1_1  // components of the objective function
  number L1_2
  number L1_3
  number L1_4  
  number L2_1
  number L2_2
  number L2_3
  number L2_4
  number L3_1
  number L3_2
  number L3_3
  number L3_4
  
  //number rc_cps_spl1
  //number rc_cps_spl2
  //number rc_cps_spl3
  //number rc_cps_spl4
  vector tacAge1(fage,lage) //TAC at age
  vector tacAge2(fage,lage)
  vector tacAge3(fage,lage)
  vector tacAge4(fage,lage)
  vector Ftac(fage,lage) //fishing mort for tac 
  vector Mvec(fage,lage)
  vector Zvec(fage,lage)
  vector Svec(fage,lage)
  
  objective_function_value negLL;


RUNTIME_SECTION
  maximum_function_evaluations 50000

PRELIMINARY_CALCS_SECTION 
  random_number_generator rnd(seed);

  //creating the natural mortality matrix; M not estimated as part of the model
  M=nat_mort;
   
  rc_cch1=rowsum(rc_catch1);                     //Total catch in a year
  rc_cch2=rowsum(rc_catch2);
  rc_cch3=rowsum(rc_catch3);
  rc_cch4=rowsum(rc_catch4);
  for(i=fyear;i<=lyear;i++)
  {
    for(j=fage;j<=lage;j++)
    {
      rc_cps1(i,j)=rc_catch1(i,j)/rc_cch1(i);     //age comp by division
      rc_cps2(i,j)=rc_catch2(i,j)/rc_cch2(i); 
      rc_cps3(i,j)=rc_catch3(i,j)/rc_cch3(i);
      rc_cps4(i,j)=rc_catch4(i,j)/rc_cch4(i);
    }   
  }

  //cout<<"Total rec catch"<<endl;   
  //cout<<rc_cch<<endl;
  //cout<<"Age composition"<<endl;
  //cout<<rc_cps<<endl;
  //exit(0);

 
  //following add by liu 
  ivector fish1(1,effSampleSize); //resampling fish from readin catch_at_age
  ivector fish2(1,effSampleSize);
  ivector fish3(1,effSampleSize);
  ivector fish4(1,effSampleSize);
  for(i=fyear;i<=lyear;i++){
    fish1.fill_multinomial(rnd,rc_cps1(i)); //based on readin age composition and effect sample size 
    fish2.fill_multinomial(rnd,rc_cps2(i)); 
    fish3.fill_multinomial(rnd,rc_cps3(i)); 
    fish4.fill_multinomial(rnd,rc_cps4(i)); 
    obsResampleCatch1(i)= CountFreqAges(fish1,fage,lage);  //total is effSampleSize
    obsResampleCatch2(i)= CountFreqAges(fish2,fage,lage);
    obsResampleCatch3(i)= CountFreqAges(fish3,fage,lage);
    obsResampleCatch4(i)= CountFreqAges(fish4,fage,lage);
    new_rc_cps1(i)=obsResampleCatch1(i)/double(effSampleSize);  //age composition based on resampling from effective sample size
    new_rc_cps2(i)=obsResampleCatch2(i)/double(effSampleSize);
    new_rc_cps3(i)=obsResampleCatch3(i)/double(effSampleSize);
    new_rc_cps4(i)=obsResampleCatch4(i)/double(effSampleSize);
    //add observation error for read in catch
    new_rc_cch1(i)=rc_cch1(i)*mfexp((stdObsCatch*randn(rnd))-(0.5*square(stdObsCatch)));  //add  lognormal err  
    new_rc_cch2(i)=rc_cch2(i)*mfexp((stdObsCatch*randn(rnd))-(0.5*square(stdObsCatch)));
    new_rc_cch3(i)=rc_cch3(i)*mfexp((stdObsCatch*randn(rnd))-(0.5*square(stdObsCatch)));
    new_rc_cch4(i)=rc_cch3(i)*mfexp((stdObsCatch*randn(rnd))-(0.5*square(stdObsCatch))); 
    new_rc_catch1(i)=new_rc_cch1(i)*new_rc_cps1(i);  //based on new obs. annual catch and new age compos., reallocate fish
    new_rc_catch2(i)=new_rc_cch2(i)*new_rc_cps2(i);
    new_rc_catch3(i)=new_rc_cch3(i)*new_rc_cps3(i);
    new_rc_catch4(i)=new_rc_cch4(i)*new_rc_cps4(i); 
  }
  //cout<<endl<<rc_cps <<endl<<endl<<new_rc_cps <<endl<<endl<<new_rc_cch<<endl;  //exit(9);
  //cout<<endl<<obsResampleCatch<<endl<<endl<<rc_catch <<endl<<endl<<new_rc_catch<<endl; exit(9);

PROCEDURE_SECTION
  get_mortality();
  get_population();
  get_catch();
  get_diffs();
  get_objective();




FUNCTION get_mortality
  qrc=mfexp(log_qrc);
  effort_devs1=mfexp(log_effort_devs1);
  effort_devs2=mfexp(log_effort_devs2);
  effort_devs3=mfexp(log_effort_devs3);
  effort_devs4=mfexp(log_effort_devs4);
  //Recreational selectivity modeled as a gamma function
  //pop1-4
  rc_s_p1=mfexp(log_rc_s_p1);
  for (j=fage; j<=lage; j++)
  {
    rc_s_age(j)=pow(j,rc_s_p1)*mfexp(-rc_s_p2*(j));  
  }

  for (i=fyear;i<=lyear;i++)
  {
  rc_s(i)=rc_s_age/rc_s_age(10);   // scaling
  }
  
  //Fishing mortality set proportional to effort
  //pop1-4
  for(i=fyear;i<=lyear;i++)
   {
    for(j=fage;j<=lage;j++)
     {
     F1(i,j)=qrc*rc_eff1(i)*rc_s(i,j)*effort_devs1(i);
     F2(i,j)=qrc*rc_eff2(i)*rc_s(i,j)*effort_devs2(i);
     F3(i,j)=qrc*rc_eff3(i)*rc_s(i,j)*effort_devs3(i);
     F4(i,j)=qrc*rc_eff4(i)*rc_s(i,j)*effort_devs4(i);
     }
   }
  Z1= M + F1;
  Z2= M + F2;
  Z3= M + F3;
  Z4= M + F4;
  S1 = mfexp(-1.0*Z1);     //Annual survival
  S2 = mfexp(-1.0*Z2); 
  S3 = mfexp(-1.0*Z3); 
  S4 = mfexp(-1.0*Z4); 

FUNCTION get_population
  //pop1 
  for (j=fage+1;j<=lage;j++)
  {
   initial1(j)=mfexp(log_initial1+initial_devs1(j));
   initial2(j)=mfexp(log_initial2+initial_devs2(j));
   initial3(j)=mfexp(log_initial3+initial_devs3(j));
   initial4(j)=mfexp(log_initial4+initial_devs4(j));
  }

  for (i=fyear;i<=lyear;i++)
  {
  n1(i,fage)=mfexp(log_mean_recruit1+recruit_devs1(i));  // yearling numbers over years
  n2(i,fage)=mfexp(log_mean_recruit2+recruit_devs2(i));
  n3(i,fage)=mfexp(log_mean_recruit3+recruit_devs3(i));
  n4(i,fage)=mfexp(log_mean_recruit4+recruit_devs4(i));
  }

  n1(fyear)(fage+1,lage)=initial1; //set abundance in the fist year for age 2 and older fish
  n2(fyear)(fage+1,lage)=initial2;
  n3(fyear)(fage+1,lage)=initial3;
  n4(fyear)(fage+1,lage)=initial4;
  
  for (i=fyear+1; i<=lyear; i++)
   {
     for (j=fage+1; j<=lage; j++)
     {  //Pop 1 mixing to area2-4 or stay in area1 are suffering the Survival rate in their destination area
     n1(i,j) =n1(i-1,j-1)*Mix(1,1)*S1(i-1,j-1)+n1(i-1,j-1)*Mix(1,2)*S2(i-1,j-1)+n1(i-1,j-1)*Mix(1,3)*S3(i-1,j-1)+n1(i-1,j-1)*Mix(1,4)*S4(i-1,j-1);
     n2(i,j) =n2(i-1,j-1)*Mix(2,1)*S1(i-1,j-1)+n2(i-1,j-1)*Mix(2,2)*S2(i-1,j-1)+n2(i-1,j-1)*Mix(2,3)*S3(i-1,j-1)+n2(i-1,j-1)*Mix(2,4)*S4(i-1,j-1);
     n3(i,j) =n3(i-1,j-1)*Mix(3,1)*S1(i-1,j-1)+n3(i-1,j-1)*Mix(3,2)*S2(i-1,j-1)+n3(i-1,j-1)*Mix(3,3)*S3(i-1,j-1)+n3(i-1,j-1)*Mix(3,4)*S4(i-1,j-1);
     n4(i,j) =n4(i-1,j-1)*Mix(4,1)*S1(i-1,j-1)+n4(i-1,j-1)*Mix(4,2)*S2(i-1,j-1)+n4(i-1,j-1)*Mix(4,3)*S3(i-1,j-1)+n4(i-1,j-1)*Mix(4,4)*S4(i-1,j-1);
     }
     //last age fish of Pop1 if survivaled in their movement destination area. will still go ack and add into last age class
     n1(i,lage)+=n1(i-1,lage)*Mix(1,1)*S1(i-1,lage)+n1(i-1,lage)*Mix(1,2)*S2(i-1,lage)+n1(i-1,lage)*Mix(1,3)*S3(i-1,lage)+n1(i-1,lage)*Mix(1,4)*S4(i-1,lage);  // the last age group is a plus age group
     n2(i,lage)+=n2(i-1,lage)*Mix(2,1)*S1(i-1,lage)+n2(i-1,lage)*Mix(2,2)*S2(i-1,lage)+n2(i-1,lage)*Mix(2,3)*S3(i-1,lage)+n2(i-1,lage)*Mix(2,4)*S4(i-1,lage); 
     n3(i,lage)+=n3(i-1,lage)*Mix(3,1)*S1(i-1,lage)+n3(i-1,lage)*Mix(3,2)*S2(i-1,lage)+n3(i-1,lage)*Mix(3,3)*S3(i-1,lage)+n3(i-1,lage)*Mix(3,4)*S4(i-1,lage);
     n4(i,lage)+=n4(i-1,lage)*Mix(4,1)*S1(i-1,lage)+n4(i-1,lage)*Mix(4,2)*S2(i-1,lage)+n4(i-1,lage)*Mix(4,3)*S3(i-1,lage)+n4(i-1,lage)*Mix(4,4)*S4(i-1,lage); 
   }
   
  for (i=fyear; i<=lyear; i++)
   {
    for (j=fage; j<=lage; j++) 
     { //abundance * weight which is proportional to SSB
     ssb1(i,j) =  n1(i,j) * wt_age(j)*mat_age(j)*propFemale;
     ssb2(i,j) =  n2(i,j) * wt_age(j)*mat_age(j)*propFemale;
     ssb3(i,j) =  n3(i,j) * wt_age(j)*mat_age(j)*propFemale;
     ssb4(i,j) =  n4(i,j) * wt_age(j)*mat_age(j)*propFemale;
     }
   } 
  
  for (i=fyear; i<=lyear; i++)
   {
    N1(i)=sum(n1(i));          // total population over years
    N2(i)=sum(n2(i));
    N3(i)=sum(n3(i));
    N4(i)=sum(n4(i));
    Total_ssb1(i) =sum(ssb1(i));
    Total_ssb2(i) =sum(ssb2(i));
    Total_ssb3(i) =sum(ssb3(i));
    Total_ssb4(i) =sum(ssb4(i));
   }

FUNCTION get_catch
//pop1-4
  n_rc1=elem_prod(elem_div(F1,Z1),elem_prod(1.0-S1,(n1*Mix(1,1)+n2*Mix(2,1)+n3*Mix(3,1)+n4*Mix(4,1))));	//Baranov catch equation for recreational harvest
  n_rc2=elem_prod(elem_div(F2,Z2),elem_prod(1.0-S2,(n1*Mix(1,2)+n2*Mix(2,2)+n3*Mix(3,2)+n4*Mix(4,2))));
  n_rc3=elem_prod(elem_div(F3,Z3),elem_prod(1.0-S3,(n1*Mix(1,3)+n2*Mix(2,3)+n3*Mix(3,3)+n4*Mix(4,3))));
  n_rc4=elem_prod(elem_div(F4,Z4),elem_prod(1.0-S4,(n1*Mix(1,4)+n2*Mix(2,4)+n3*Mix(3,4)+n4*Mix(4,4))));
  //d_M1 =elem_prod(elem_div(M,Z1),elem_prod(1.0-S1,(n1*Mix(1,1)+n2*Mix(2,1)+n3*Mix(3,1)+n4*Mix(4,1))));     //Deaths due to nat mort(M)   , not used anymore
  //d_M2 =elem_prod(elem_div(M,Z2),elem_prod(1.0-S2,(n1*Mix(1,2)+n2*Mix(2,2)+n3*Mix(3,2)+n4*Mix(4,2))));
  //d_M3 =elem_prod(elem_div(M,Z3),elem_prod(1.0-S3,(n1*Mix(1,3)+n2*Mix(2,3)+n3*Mix(3,3)+n4*Mix(4,3)))); 
  //d_M4 =elem_prod(elem_div(M,Z4),elem_prod(1.0-S4,(n1*Mix(1,4)+n2*Mix(2,4)+n3*Mix(3,4)+n4*Mix(4,4))));
  //Estimates total catch by year
  for (i=fyear; i<=lyear; i++)
   {
   N_rc1(i)=sum(n_rc1(i));
   N_rc2(i)=sum(n_rc2(i));
   N_rc3(i)=sum(n_rc3(i));
   N_rc4(i)=sum(n_rc4(i));   
   }

   //Estimated age composition of recreational catch
  for (i=fyear; i<=lyear; i++)
  {
    n_rc_cps1(i)=n_rc1(i)/N_rc1(i);
    n_rc_cps2(i)=n_rc2(i)/N_rc2(i);
    n_rc_cps3(i)=n_rc3(i)/N_rc3(i);
    n_rc_cps4(i)=n_rc4(i)/N_rc4(i);
  }


FUNCTION get_diffs
//calculate differences in total catch, survvey cpe, and juvenile index
  for(i=fyear;i<=lyear;i++)
  {
      rc_diff1(i)=(log(new_rc_cch1(i)+0.000001)-log(N_rc1(i)+0.0000001)); 
      rc_diff2(i)=(log(new_rc_cch2(i)+0.000001)-log(N_rc2(i)+0.0000001)); 
      rc_diff3(i)=(log(new_rc_cch3(i)+0.000001)-log(N_rc3(i)+0.0000001)); 
      rc_diff4(i)=(log(new_rc_cch4(i)+0.000001)-log(N_rc4(i)+0.0000001)); 
      re1(i)=(N_rc1(i)+0.000001-new_rc_cch1(i)+0.000001)/(new_rc_cch1(i)+0.000001);
      re2(i)=(N_rc2(i)+0.000001-new_rc_cch2(i)+0.000001)/(new_rc_cch2(i)+0.000001);
      re3(i)=(N_rc3(i)+0.000001-new_rc_cch3(i)+0.000001)/(new_rc_cch3(i)+0.000001);
      re4(i)=(N_rc4(i)+0.000001-new_rc_cch4(i)+0.000001)/(new_rc_cch4(i)+0.000001);
      abs_re1(i)=re1(i);
      abs_re2(i)=re2(i);
      abs_re3(i)=re3(i);
      abs_re4(i)=re4(i);
      if(re1(i)<0) abs_re1(i)*=-1;
      if(re2(i)<0) abs_re2(i)*=-1;
      if(re3(i)<0) abs_re3(i)*=-1;
      if(re4(i)<0) abs_re4(i)*=-1;
 }
 mre1=quantile(re1,0.5);
 mre2=quantile(re2,0.5);
 mre3=quantile(re3,0.5);
 mre4=quantile(re4,0.5);
 mare1=quantile(abs_re1,0.5);
 mare2=quantile(abs_re2,0.5);
 mare3=quantile(abs_re3,0.5);
 mare4=quantile(abs_re4,0.5);
 
FUNCTION get_objective
  variance_ratio=0.25;
  log_sigma_effort=log(sqrt((1./variance_ratio)*square(mfexp(log_sigma_rc))));
  L1_1 = size_count(rc_diff1)*log_sigma_rc+1.0/(2.0*square(mfexp(log_sigma_rc)))*norm2(rc_diff1) ;       //Likelihood for recreati 
  L1_2 = size_count(rc_diff2)*log_sigma_rc+1.0/(2.0*square(mfexp(log_sigma_rc)))*norm2(rc_diff2) ;
  L1_3 = size_count(rc_diff3)*log_sigma_rc+1.0/(2.0*square(mfexp(log_sigma_rc)))*norm2(rc_diff3) ;
  L1_4 = size_count(rc_diff4)*log_sigma_rc+1.0/(2.0*square(mfexp(log_sigma_rc)))*norm2(rc_diff4) ;   
  L2_1 = -sum(double(effSampleSize)*elem_prod(rc_cps1,log(0.0000001+n_rc_cps1)));     //age comp for recreational fishery, old way
  L2_2 = -sum(double(effSampleSize)*elem_prod(rc_cps2,log(0.0000001+n_rc_cps2)));
  L2_3 = -sum(double(effSampleSize)*elem_prod(rc_cps3,log(0.0000001+n_rc_cps3)));
  L2_4 = -sum(double(effSampleSize)*elem_prod(rc_cps4,log(0.0000001+n_rc_cps4)));
  L3_1 = size_count(log_effort_devs1)*log_sigma_effort+1.0/(2.0*square(mfexp(log_sigma_effort)))*norm2(log_effort_devs1);
  L3_2 = size_count(log_effort_devs2)*log_sigma_effort+1.0/(2.0*square(mfexp(log_sigma_effort)))*norm2(log_effort_devs2);
  L3_3 = size_count(log_effort_devs3)*log_sigma_effort+1.0/(2.0*square(mfexp(log_sigma_effort)))*norm2(log_effort_devs3);
  L3_4 = size_count(log_effort_devs4)*log_sigma_effort+1.0/(2.0*square(mfexp(log_sigma_effort)))*norm2(log_effort_devs4);
 
  negLL=L1_1 +L1_2 +L1_3 +L1_4 +L2_1 +L2_2 +L2_3 +L2_4 +L3_1 +L3_2 +L3_3 +L3_4;
  
  
//only being called in report section
FUNCTION get_TAC
  dvar_matrix abund1(lyear,lyear+2,fage,lage); //only used for tac calculation
  dvar_matrix abund2(lyear,lyear+2,fage,lage); 
  dvar_matrix abund3(lyear,lyear+2,fage,lage); 
  dvar_matrix abund4(lyear,lyear+2,fage,lage); 
  
  abund1(lyear)=n1(lyear); //fill in the last yr abundance from n
  abund2(lyear)=n2(lyear);
  abund3(lyear)=n3(lyear);
  abund4(lyear)=n4(lyear);
  
  //cout<<abund(lyear)<<endl;

  
  //use last-year Z as the Z for year last+1  ,use average recruitment
  for(j=fage+1; j<=lage; j++) {
   abund1(lyear+1,j) =abund1(lyear,j-1)*Mix(1,1)*mfexp(-Z1(lyear,j-1))+ abund1(lyear,j-1)*Mix(1,2)*mfexp(-Z2(lyear,j-1))+ abund1(lyear,j-1)*Mix(1,3)*mfexp(-Z3(lyear,j-1))+ abund1(lyear,j-1)*Mix(1,4)*mfexp(-Z4(lyear,j-1));  
   abund2(lyear+1,j) =abund2(lyear,j-1)*Mix(2,1)*mfexp(-Z1(lyear,j-1))+ abund2(lyear,j-1)*Mix(2,2)*mfexp(-Z2(lyear,j-1))+ abund2(lyear,j-1)*Mix(2,3)*mfexp(-Z3(lyear,j-1))+ abund2(lyear,j-1)*Mix(2,4)*mfexp(-Z4(lyear,j-1));
   abund3(lyear+1,j) =abund3(lyear,j-1)*Mix(3,1)*mfexp(-Z1(lyear,j-1))+ abund3(lyear,j-1)*Mix(3,2)*mfexp(-Z2(lyear,j-1))+ abund3(lyear,j-1)*Mix(3,3)*mfexp(-Z3(lyear,j-1))+ abund3(lyear,j-1)*Mix(3,4)*mfexp(-Z4(lyear,j-1));
   abund4(lyear+1,j) =abund4(lyear,j-1)*Mix(4,1)*mfexp(-Z1(lyear,j-1))+ abund4(lyear,j-1)*Mix(4,2)*mfexp(-Z2(lyear,j-1))+ abund4(lyear,j-1)*Mix(4,3)*mfexp(-Z3(lyear,j-1))+ abund4(lyear,j-1)*Mix(4,4)*mfexp(-Z4(lyear,j-1));      
  }
  
  abund1(lyear+1,lage)+=abund1(lyear,lage)*Mix(1,1)*mfexp(-Z1(lyear,lage))+ abund1(lyear,lage)*Mix(1,2)*mfexp(-Z2(lyear,lage))+ abund1(lyear,lage)*Mix(1,3)*mfexp(-Z3(lyear,lage))+ abund1(lyear,lage)*Mix(1,4)*mfexp(-Z4(lyear,lage));  // the last age group is a plus age group
  abund2(lyear+1,lage)+=abund2(lyear,lage)*Mix(2,1)*mfexp(-Z1(lyear,lage))+ abund2(lyear,lage)*Mix(2,2)*mfexp(-Z2(lyear,lage))+ abund2(lyear,lage)*Mix(2,3)*mfexp(-Z3(lyear,lage))+ abund2(lyear,lage)*Mix(2,4)*mfexp(-Z4(lyear,lage));
  abund3(lyear+1,lage)+=abund3(lyear,lage)*Mix(3,1)*mfexp(-Z1(lyear,lage))+ abund3(lyear,lage)*Mix(3,2)*mfexp(-Z2(lyear,lage))+ abund3(lyear,lage)*Mix(3,3)*mfexp(-Z3(lyear,lage))+ abund3(lyear,lage)*Mix(3,4)*mfexp(-Z4(lyear,lage));
  abund4(lyear+1,lage)+=abund4(lyear,lage)*Mix(4,1)*mfexp(-Z1(lyear,lage))+ abund4(lyear,lage)*Mix(4,2)*mfexp(-Z2(lyear,lage))+ abund4(lyear,lage)*Mix(4,3)*mfexp(-Z3(lyear,lage))+ abund4(lyear,lage)*Mix(4,4)*mfexp(-Z4(lyear,lage));
   
  abund1(lyear+1,fage)=sum(column(n1,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);    //use average recruitment
  abund2(lyear+1,fage)=sum(column(n2,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);    //use average recruitment
  abund3(lyear+1,fage)=sum(column(n3,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);    //use average recruitment
  abund4(lyear+1,fage)=sum(column(n4,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);    //use average recruitment
  
  dvariable tmpMort1;
  dvariable tmpMort2;
  dvariable tmpMort3;
  dvariable tmpMort4; 
  //use average Z as the Z for year last+1  ,use average recruitment
  for(j=fage+1; j<=lage; j++){
   tmpMort1=mean(column(Z1,j-1)(lyear-lastNumYrMort+1,lyear) ); //how many last yrs mortality as mean mortality for TAC
   tmpMort2=mean(column(Z2,j-1)(lyear-lastNumYrMort+1,lyear) ); //how many last yrs mortality as mean mortality for TAC
   tmpMort3=mean(column(Z3,j-1)(lyear-lastNumYrMort+1,lyear) ); //how many last yrs mortality as mean mortality for TAC
   tmpMort4=mean(column(Z4,j-1)(lyear-lastNumYrMort+1,lyear) ); //how many last yrs mortality as mean mortality for TAC   
   abund1(lyear+2,j) =abund1(lyear+1,j-1)*Mix(1,1)*mfexp(-tmpMort1)+ abund1(lyear+1,j-1)*Mix(1,2)*mfexp(-tmpMort2)+ abund1(lyear+1,j-1)*Mix(1,3)*mfexp(-tmpMort3)+ abund1(lyear+1,j-1)*Mix(1,4)*mfexp(-tmpMort4);  
   abund2(lyear+2,j) =abund2(lyear+1,j-1)*Mix(2,1)*mfexp(-tmpMort1)+ abund2(lyear+1,j-1)*Mix(2,2)*mfexp(-tmpMort2)+ abund2(lyear+1,j-1)*Mix(2,3)*mfexp(-tmpMort3)+ abund2(lyear+1,j-1)*Mix(2,4)*mfexp(-tmpMort4);
   abund3(lyear+2,j) =abund3(lyear+1,j-1)*Mix(3,1)*mfexp(-tmpMort1)+ abund3(lyear+1,j-1)*Mix(3,2)*mfexp(-tmpMort2)+ abund3(lyear+1,j-1)*Mix(3,3)*mfexp(-tmpMort3)+ abund3(lyear+1,j-1)*Mix(3,4)*mfexp(-tmpMort4);
   abund4(lyear+2,j) =abund4(lyear+1,j-1)*Mix(4,1)*mfexp(-tmpMort1)+ abund4(lyear+1,j-1)*Mix(4,2)*mfexp(-tmpMort2)+ abund4(lyear+1,j-1)*Mix(4,3)*mfexp(-tmpMort3)+ abund4(lyear+1,j-1)*Mix(4,4)*mfexp(-tmpMort4);         
  }
  tmpMort1=mean(column(Z1,lage)(lyear-lastNumYrMort+1,lyear) );
  tmpMort2=mean(column(Z2,lage)(lyear-lastNumYrMort+1,lyear) );
  tmpMort3=mean(column(Z3,lage)(lyear-lastNumYrMort+1,lyear) );
  tmpMort4=mean(column(Z4,lage)(lyear-lastNumYrMort+1,lyear) );
  
  abund1(lyear+2,lage)+=abund1(lyear+1,lage)*Mix(1,1)*mfexp(-tmpMort1)+ abund1(lyear+1,lage)*Mix(1,2)*mfexp(-tmpMort2)+ abund1(lyear+1,lage)*Mix(1,3)*mfexp(-tmpMort3)+ abund1(lyear+1,lage)*Mix(1,4)*mfexp(-tmpMort4);  // the last age group is a plus age group
  abund2(lyear+2,lage)+=abund2(lyear+1,lage)*Mix(2,1)*mfexp(-tmpMort1)+ abund2(lyear+1,lage)*Mix(2,2)*mfexp(-tmpMort2)+ abund2(lyear+1,lage)*Mix(2,3)*mfexp(-tmpMort3)+ abund2(lyear+1,lage)*Mix(2,4)*mfexp(-tmpMort4);
  abund3(lyear+2,lage)+=abund3(lyear+1,lage)*Mix(3,1)*mfexp(-tmpMort1)+ abund3(lyear+1,lage)*Mix(3,2)*mfexp(-tmpMort2)+ abund3(lyear+1,lage)*Mix(3,3)*mfexp(-tmpMort3)+ abund3(lyear+1,lage)*Mix(3,4)*mfexp(-tmpMort4);
  abund4(lyear+2,lage)+=abund4(lyear+1,lage)*Mix(4,1)*mfexp(-tmpMort1)+ abund4(lyear+1,lage)*Mix(4,2)*mfexp(-tmpMort2)+ abund4(lyear+1,lage)*Mix(4,3)*mfexp(-tmpMort3)+ abund4(lyear+1,lage)*Mix(4,4)*mfexp(-tmpMort4);
  
  abund1(lyear+2,fage)=sum(column(n1,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec);  //the last 5yrs first age as mean recruitment
  abund2(lyear+2,fage)=sum(column(n2,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec); 
  abund3(lyear+2,fage)=sum(column(n3,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec); 
  abund4(lyear+2,fage)=sum(column(n4,fage)(lyear-lastNumYrRec+1,lyear) )/double(lastNumYrRec); 
  
  Ftac=(tot_targ_Z-nat_mort)*rc_s(lyear);  //use target total mortality and last yr selectivity
  Mvec=M(lyear);
  Zvec=Mvec+Ftac;
  Svec=mfexp(-1.0*Zvec);
  
  tacAge1=elem_prod(elem_div(Ftac,Zvec),elem_prod(1.0-Svec,(abund1(lyear+2)*Mix(1,1)+abund2(lyear+2)*Mix(2,1)+abund3(lyear+2)*Mix(3,1)+abund4(lyear+2)*Mix(4,1))));
  tacAge2=elem_prod(elem_div(Ftac,Zvec),elem_prod(1.0-Svec,(abund1(lyear+2)*Mix(1,2)+abund2(lyear+2)*Mix(2,2)+abund3(lyear+2)*Mix(3,2)+abund4(lyear+2)*Mix(4,2))));
  tacAge3=elem_prod(elem_div(Ftac,Zvec),elem_prod(1.0-Svec,(abund1(lyear+2)*Mix(1,3)+abund2(lyear+2)*Mix(2,3)+abund3(lyear+2)*Mix(3,3)+abund4(lyear+2)*Mix(4,3))));
  tacAge4=elem_prod(elem_div(Ftac,Zvec),elem_prod(1.0-Svec,(abund1(lyear+2)*Mix(1,4)+abund2(lyear+2)*Mix(2,4)+abund3(lyear+2)*Mix(3,4)+abund4(lyear+2)*Mix(4,4))));
  
  
//return num of fish at same age 
FUNCTION dvector CountFreqAges(const ivector& age,int fage,int lage)
  int numAge=lage-fage+1;
  ivector uniqAge(fage,lage);
  uniqAge.fill_seqadd(fage,1); 
  dvector freq(fage,lage); 
  freq.initialize();    

  //count the number of ages which at the same age
  for(int i=age.indexmin();i<=age.indexmax();i++){
    for(int j=fage;j<=lage;j++){
      if(age(i)==uniqAge(j)){ // found match
        freq(j)+=1;
        break;
      }
    }
  }
  return freq;



REPORT_SECTION
  get_TAC();
  report<<sum(tacAge1)<<endl;
  report<<sum(tacAge2)<<endl;
  report<<sum(tacAge3)<<endl;
  report<<sum(tacAge4)<<endl;
  report<<objective_function_value::gmax <<endl;
  report<<N1(lyear)<<endl;  //estimated total abundance in the last year of the assessment  
  report<<N2(lyear)<<endl;
  report<<N3(lyear)<<endl;
  report<<N4(lyear)<<endl;
  report<<mre1<<endl;  //median relative error between observed and predicted catch for the 20 yrs included in this assessment
  report<<mre2<<endl;
  report<<mre3<<endl;
  report<<mre4<<endl;
  report<<mare1<<endl;  //median absolute relative error between observed and predicted for the 20yrs included in this assessment
  report<<mare2<<endl;
  report<<mare3<<endl;
  report<<mare4<<endl;
  report<<Total_ssb1(lyear)<<endl;
  report<<Total_ssb2(lyear)<<endl;
  report<<Total_ssb3(lyear)<<endl;
  report<<Total_ssb4(lyear)<<endl;
  /*
  report << "Relative errir" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << re(i) << endl;
  }
  report << "Absolute Relative errir" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << abs_re(i) << endl;
  }
  report << "Total Abundance" << endl;
  for (i=fyear;i<=lyear;i++)
  {
  report << N(i) << endl;
  }
  report << "Abundance at Age" << endl;
  report << n << endl;
  report << "Fishing Mortality" << endl;
  report << F << endl;
  report << "Natural Mortality" << endl;
  report << M << endl;
  report << "Total Mortality" << endl;
  report << Z << endl;
  report << "Observed Recreational Catch" << " " << "Predicted Recreational Catch" << endl;
  for (i=fyear;i<=lyear;i++)
      {
          report << new_rc_cch(i) << "                 " << N_rc(i) << endl; //change by liu
      }

  report << "Recreational Selectivity" << endl;
  report << rc_s<< endl;  //changed by liu
  report << "Recreational catch standard deviation" << endl;  
  report << mfexp(log_sigma_rc) << endl;
  report << "Recreational catchability" << endl;
  report << qrc << endl;
  report << "Recreational differences" << endl;
  report << rc_diff << endl;
  report << "F_TAC at Age" << endl;
  report << Ftac << endl;
  report << "TAC at Age" << endl;
  report << tacAge << endl;
  */





