/* Work in progress*/
/*******************************************************/
/* Proportional hazards regression in SAS              */
/*******************************************************/
/* The following is code I wrote for a presentation in */
/* Survival Analysis. It contains a proportional       */
/* hazards regression analysis carried out in SAS.     */
/*                                                     */
/* The dataset contains information about 137 lung     */
/* cancer patients with four disparate types of cancer */ 
/* cells. Each patient was treated with a standard     */
/* therapy or an experimental one. Several covariates  */
/* are included.                                       */   
/*******************************************************/










/*******************************************************/
/*Import data and preprocessing steps                  */
/*******************************************************/


filename csvFile 
	url "https://github.com/ypaulsen/Survival-Analysis/raw/main/valung.csv";

proc import datafile=csvFile 
	out=valung replace dbms=csv; 
run;

/*Data structure:*/ 
/* 
treatment: therapy 
cancer type: cell = {Squamous, Small, Large, Adeno} 
time: t {time in days}
outcome: dead {dead, censored} 
Covariates: 
  Numeric: 
    kps diagtime age 
  Other: 
    prior {Yes No}
*/

/*Look at data*/
proc print data=valung; 
run; 	

/*Rename*/ 
data lung; 
	set valung; 
run; 

/*Code categorical variables as integers*/ 

data lung;
	set lung; 
	dead_int  = input(dead_int, 1.);   /*Create new variables*/
	prior_int  = input(prior_int, 1.);
	therapy_int  = input(therapy_int, 1.);
	cell_int  = input(cell_int, 1.); 
	if dead = 'dead' then dead_int = 1;      /*Assign values to new variables*/
	else dead_int = 0;
	if prior = 'yes' then prior_int = 1;      
	else prior_int = 0;

	if therapy = 'standard' then therapy_int = 0; 
	else therapy_int=1;
	if cell = 'Squamous' then cell_int = 1; 
	else if cell = 'Small' then cell_int = 2; 
	else if cell = 'Adeno' then cell_int = 3; 
	else cell_int = 4;
run;

/*Look at data*/
proc print data=lung; 
run;

/*******************************************************/
/*  End of data processing steps                       */  
/*******************************************************/







/*******************************************************/
/* PROC Lifetest                                       */ 
/*******************************************************/
/* Preliminary analyses                                */
/*******************************************************/



/* Data:  
treatment: therapy 
z: time: t 
delta: outcome: dead_int 
cancer type: cell
Coviariates: kps diagtime age prior
*/




/*******************/
/* testing therapy */
/*******************/

/* Kaplan Meier survivor estimates with                */
/* Nelson Aalen cummulative hazard functions           */
proc lifetest data=lung method=km nelson plots=(survival(cl), ls, lls) 
	outsurv=a1; 
/* lls for proportional hazards */
/* ls for cummulative hazard */
	time t*dead_int(0);
	strata therapy; 
run;




























/* Testing all covariates with forward stepwise        */
/* elimination.                                        */
proc lifetest data=lung method=km plots = (survival(cl), ls, lls); 
	time t*dead_int(0);
	test kps diagtime age; 
run; 

/*Only kps is significant???                           */



proc lifetest data=lung method=km plots = (survival(cl), ls, lls); 
	time t*dead_int(0);
	strata cell; 
run; 






























/* Kaplan Meier survivor estimates linear confidence   */
proc lifetest data=lung method=km nelson conftype=linear plots=(survival(cl), ls, lls) 
	outsurv=a1; 
	time t*dead_int(0);
	strata therapy; 
run;



/*******************************************************/
/* testing therapy within cell types                   */
/*******************************************************/


/* Kaplan Meier survivor estimates*/   
proc lifetest data=lung method=km nelson plots=(survival(cl), ls, lls)
	outsurv=a2; 
	time t*dead_int(0);
	strata cell/ group = therapy; 
run;


/*******************************************************/
/* End of proc lifetest section.                       */
/*******************************************************/





/*******************************************************/
/* Graphically checking the distribution of Y.         */
/*                                                     */
/*******************************************************/      
/* Most linear plot == best model for distribution     */
/* The code is included here for example puroses.      */
/*******************************************************/

/*Generate data*/  
data a3;
	set a1;
	s = survival;
	logH = log(-log(s));
	lnorm = probit(1-s);
	logit = log(s/(1-s));
	ltime = log(t);
run;


/*logit for log-logistic, logH for weibull and lnorm for log-normal distribution */
proc gplot data=a3;
	symbol1 i=join width=2 value=triangle c=steelblue;
	symbol2 i=join width=2 value=circle c=grey;
	plot logit*ltime=kps logH*ltime=kps lnorm*ltime=kps; 
run;



proc gplot data=a3;
	symbol1 i=join width=2 value=triangle c=steelblue;
	symbol2 i=join width=2 value=circle c=grey;
	plot logit*ltime=therapy logH*ltime=therapy lnorm*ltime=therapy; 
run;


proc gplot data=a3;
title "Graphically checking for proportional hazards property";
plot logH*t=therapy;
run;

/*******************************************************/
/* End of cum hazard plots                             */
/*******************************************************/






/*******************************************************/
/*PROC PHREG                                           */
/*******************************************************/

/* Full model with 2 methods for partial likelihood    */
/* estimation.                                         */
/* ties = Breslow is default, but generally inferior   */
/* for heavily tied data.                              */
/* Fit with all three: if similar, then the data are   */
/* not heavily tied.                                   */ 
/* If it fails to converge (not the case here) then    */
/* look at preliminary analyses to find out which      */
/* features are unimportant -> drop those and it may   */
/* converge.                                           */ 
title;
proc phreg data=lung;
	class therapy_int cell_int;
    model t*dead_int(0) = therapy_int kps diagtime age prior_int cell_int/ 
	ties=efron;
run;


proc phreg data=lung;
	class cell therapy;
    model t*dead_int(0) = therapy kps diagtime age prior_int cell/ 
	ties=exact;
run;

/* Full model, efron method, with backwards selection. */
/* Conduct model selection with efron method to save   */
/* computation time on large datasets.                 */ 

proc phreg data=lung;
	class cell therapy prior;
    model t*dead_int(0)= therapy kps diagtime age prior cell/ 
	ties=efron selection=backward;
run;



/* Backward selection eleminated 'therapy' with p=.19 */ 
/* But it is the feature of interest so I will leave  */
/* it in to test in my final model.                   */ 

/* Fit the *final* model with exact method.           */
/*(Including therapy)                                 */ 

proc phreg data=lung;
	class cell therapy;
    model t*dead_int(0)=therapy kps cell/ 
	ties=exact;
run;



/* Fit the *final* model with exact method.           */
/*(Excluding therapy)                                 */ 

proc phreg data=lung;
	class cell;
    model t*dead_int(0)= kps cell/ 
	ties=exact;
run;

 
/* Baseline survival function 
/* Finish this section later 

data null;
	input kps cell;
	cards;
0 0 
run;

proc phreg data=lung;
	model t*dead_int(0)=kps cell
	/ties=exact covb;
 	baseline out=a covariates=null survival=s lower=lcl upper=ucl
	cumhaz=H lowercumhaz=lH uppercumhaz=uH;
run;
*/
/*******************************************************/






/*******************************************************/
/* checking proportional hazard assumption             */
/* with resampling                                     */
/*******************************************************/


/* Including therapy                                   */   

proc phreg data=lung;
	class cell therapy; 
	model t*dead_int(0)= therapy kps cell/ ties=exact;
	assess ph/ resample;
run;



/* Excluding therapy                                   */   

proc phreg data=lung;
	class cell; 
	model t*dead_int(0)= kps cell/ ties=exact;
	assess ph/ resample;
run;





