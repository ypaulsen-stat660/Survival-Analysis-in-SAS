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

/*Change some variables for analysis.*/ 
/*Create new variable called 'dead_int' with "dead" recorded as 1 and "censored" as 0.*/
/*Create new variable called 'prior_int' with "Yes" recorded as 1 and "No" as 0.*/
data lung;
	set lung; 
	dead_int  = input(dead_int, 1.);   /*Create new variables*/
	prior_int  = input(prior_int, 1.); 
	if dead = 'dead' then dead_int = 1;      /*Assign values to new variables*/
	else dead_int = 0;
	if prior = 'yes' then prior_int = 1;      
	else prior_int = 0;
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


/* Data:  
treatment: therapy 
z: time: t 
delta: outcome: dead_int 
cancer type: cell
Coviariates: kps diagtime age prior
*/




/*******************/
/* testing therapy */


/* Kaplan Meier survivor estimates*/
proc lifetest data=lung method=km plots=(survival(cl), lls) outsurv=a1; 
/* lls for proportional hazards */
	time t*dead_int(0);
	strata therapy; 
run;


/* Kaplan Meier survivor estimates*/
proc lifetest data=lung method=nelson plots=ls; 
/* ls for cummulative hazard*/
	time t*dead_int(0);
	strata therapy; 
run;







/*************************************/
/* testing therapy within cell types */

/* Kaplan Meier survivor estimates*/   
proc lifetest data=lung method=km plots=(survival(cl), lls)
	outsurv=a2; 
	time t*dead_int(0);
	strata cell/ group=therapy; 
run;


/* Nelson Aalen Hazard Plots*/   
proc lifetest data=lung method=nelson plots=ls; 
	time t*dead_int(0);
	strata cell/ group=therapy; 
run;


/*******************************************************/
/* End of proc lifetest section.                       */
/*******************************************************/



/*******************************************************/
/* Graphically checking the distribution of Y.         */
/*                                                     */
/*******************************************************/
/* Below analyses are unnecessary in this case since   */
/* the plots above show that exp works here.           */
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

data a4;
	set a2;
	s = survival;
	logH = log(-log(s));
	lnorm = probit(1-s);
	logit = log(s/(1-s));
	ltime = log(t);
run;

*proc print data=a2; 
*run; 

/*logit for log-logistic, logH for weibull and lnorm for log-normal distribution */
proc gplot data=a3;
	symbol1 i=join width=2 value=triangle c=steelblue;
	symbol2 i=join width=2 value=circle c=grey;
	plot logit*ltime=therapy logH*ltime=therapy lnorm*ltime=therapy; 
run;



proc gplot data=a4;
	symbol1 i=join width=2 value=triangle c=steelblue;
	symbol2 i=join width=2 value=circle c=grey;
	plot logit*ltime=therapy logH*ltime=therapy lnorm*ltime=therapy; 
run;


/*******************************************************/
/* End of cum hazard plots                             */
/*******************************************************/






/*******************************************************/
/*PROC PHREG                                           */
/*******************************************************/

/* Full model with 2 methods for ties                  */
/* Breslow is default, but generally inferior          */

proc phreg data=lung;
	class cell therapy;
    model t*dead_int(0) = therapy kps diagtime age prior_int cell/ 
	ties=efron;
run;

proc phreg data=lung;
	class cell;
    model t*dead_int(0) = kps diagtime age prior_int cell/ 
	ties=exact;
run;

/* Full model, efron method, with backwards selection. */
/* Conduct model selection with efron method to save   */
/* computation time on large datasets.                 */ 

proc phreg data=lung;
	class cell therapy prior_int;
    model t*dead_int(0)= therapy kps diagtime age prior_int cell/ 
	ties=efron selection=backward;
run;



/* Backward selection eleminated 'therapy' with p=.19 */ 
/* But it is the feature of interest so I will leave  */
/* to test it in my final model.                      */ 

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


/*a small change */



