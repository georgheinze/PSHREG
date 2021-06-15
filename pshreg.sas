/*
SAS macro PSHREG: prepares data and fits a proportional subdistribution hazards model as proposed by Fine and Gray (1999)
Copyright (C) Georg Heinze (georg.heinze@meduniwien.ac.at)

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License (http://www.gnu.org/licenses/gpl-2.0.txt) for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

%macro pshreg(data=_last_, time=, cens=, failcode=1, cencode=0, varlist=, class=, cengroup=, firth=0, 
options=, id=, finegray=, cuminc=0, by=, censcrr=_censcrr_, out=dat_crr, weights=0, call=_pshregopt, clean=1, action=estimate,
missing=drop, statustab=1, delwork=1, tiedcens=after, admin=);

*** please send any comments to: georg.heinze@meduniwien.ac.at;
*** the macro and a User`s Guide are available for download at: 
	http://cemsiis.meduniwien.ac.at/en/kb/science-research/software/statistical-software/pshreg/ ;
*** For version >=2013.09, I appreciate some comments from Lovedeep Gondara, BC Cancer Agency, CDN
    which led to addition of the class option and sorting of input data in case of BY processing;
*** For version >=2014.06, I acknowledge the suggestions of Jerome Rapion and Akos F. Pap who discovered
    an issue when estimating a stratified Fine-Gray model.;
*** For version >=2014.09, I acknowledge the suggestions of an anonymous reviewer who suggested to account for the possibility
    of purely administrative censoring that led to the introduction of the admin option.;

%let version=2014.09;
%let build=201409301506;

/*
 data        ...	SAS data set
 time        ...	time variable
 cens        ...	censoring variable, with codes for:
                 		censoring (specified in cencode), 
                 		event of interest (specified in failcode), 
                 		competing event (all other values);
 failcode=1  ... 	code of cens for event of interest
 cencode=0   ...	code of cens for censoring
 varlist     ...	variable(s) to be used as covariates or as strata
 class       ...    defines variable(s) as factors (like CLASS statement of SAS)
 options     ...	any options passed to proc phreg (e.g., rl=pl)
 firth       ...    request Firths bias correction
 id          ...	patient identifier (not necessarily needed if there is one line per patient)
 out=dat_crr ... 	output data set, which is the input data set modified for analysis by Fine-Gray model
 by          ...	variable to define subset (for efficient processing of multiple data sets of the same structure)
 cuminc=0    ... 	to plot cumulative incidence curves by empirical subdistribution hazard method 
                    (stratified by the levels of the first variable specified in varlist);
 call        ...    name of data set with assigned values of macro options
 clean       ...	1 (default) if working data set should be cleaned,i.e., keeping only relevant variables mentioned in the macro call
 action      ...    action=estimate (default) to run PROC PHREG on the modified data set and thus compute the Fine-Gray model,
                    action=code to only compute the modified data set and put the PROC PHREG call into the Log file 
                    (replaces the legacy finegray option)
 missing     ...    missing=keep (default) to keep missing values in the modified data set
                    missing=drop to delete observations with missing values in any variables of the varlist option
 delwork     ...    1 (default) to delete all working data sets on exit
 tiedcens    ...    how censored times that are tied with event times should be handled: after event times (default), or before?
 admin       ...    administrative censoring time (if there is purely administrative censoring)
 

 *****************************************************************************************************************************************;
 The following option is experimental, and their potential in real data analyses has not yet been investigated:                    ;
 *****************************************************************************************************************************************;
 weights=1    ... 	to apply weights to the risk set according to the inverse probability of being uncensored
         2    ...    applies weights to estimate average subdistribution hazard ratios according to Schemper, Wakounig and Heinze (2009);
         
 *****************************************************************************************************************************************;
 legacy options (still available, but with no default value):
 finegray=1  ... 	to estimate the Fine-Gray proportional subdistribution hazards model using as 
					covariates all variables specified in varlist;
 
 *****************************************************************************************************************************************;
 *** A typical macro call to fit a Fine-Gray model:                                                                                *******;
     assuming data is saved in a data set called mydata,
     time_variable is the time variable,
     status_variable is the status variable, with codes: 0 for censored, 1 for event of interest, 2 for competing event;
     X1 X2 X3 are three covariates;

     %pshreg(data=mydata, time=time_variable, cens=status_variable, varlist=X1 X2 X3);
 
 *****************************************************************************************************************************************;
 */

*** finegray is a legacy option, the new version uses the more intuitive option "action";

%if &finegray=1 %then %let action=estimate;
%if &finegray=0 %then %let action=code;


* this little code chunk extracts variable names out of the variable list;

%let nvar=0;
%do %while(%scan(&varlist,&nvar+1)~=);
 %let nvar=%eval(&nvar+1);
 %let var&nvar=%scan(&varlist,&nvar);
%end;

%let nclass=0;
%if &class ne %then %do;
	%do %while(%scan(&class,&nclass+1)~=);
 	%let nclass=%eval(&nclass+1);
 	%let class&nclass=%scan(&class,&nclass);
	%end;
%end;


*** preparing a data set which holds all specifications for the macro option;

data &call;
length macrovariable $13 value $12 remark $ 49;
label macrovariable="Macro option" value="Assigned value" remark="Remark";
macrovariable="data"; value="&data"; remark="Input data set"; output;
macrovariable="time"; value="&time"; remark="Time variable"; output;
macrovariable="cens"; value="&cens"; remark="Censoring variable"; output;
macrovariable="failcode"; value="&failcode"; remark="Code for event of interest"; output;
macrovariable="cencode"; value="&cencode"; remark="Code for censored observation"; output;
macrovariable="tiedcens"; value="&tiedcens"; remark="How censored times tied with event times"; output;
macrovariable=""; value=""; remark="should be treated"; output;
macrovariable="admin"; value="&admin"; remark="Administrative censoring time variable"; output;
macrovariable="varlist";
value="&varlist";
remark="List of covariables";
%if %length(&varlist) le 12 %then %do;
	output;
	%end;
%else %do;
	%do j=1 %to &nvar;
		%if &j ne 1 %then %do;
			macrovariable=""; remark="";
		%end;
		value="&&var&j";
		output;
	%end;
%end;
macrovariable="class"; value="&class"; remark="List of class variables"; 
	%if %eval(&nclass) le 1  %then %do;
		output;
	%end;
	%else %do;
		%do j=1 %to &nclass;
			%if &j eq 1 %then %do; 
				macrovariable="class"; remark="List of class variables";
			%end;
			%else %do;
				macrovariable=""; remark="";
			%end;
			value="&&class&j";
			output;
		%end;
	%end;
macrovariable="options"; value="&options"; remark="Options to be passed to PROC PHREG"; output; remark="";
macrovariable="firth"; value="&firth"; 
	if &firth=1 then remark="Firth`s bias correction"; 
	if &firth=0 then remark="Standard ML estimation, no Firth correction"; 
	output; remark="";
if &firth=1 then do;
	macrovariable=""; value=""; remark="  (disables robust covariance estimates (COVS))"; output;
end;
macrovariable="id"; value="&id"; remark="Subject identifier"; output; 
macrovariable="by"; value="&by"; remark="BY processing variable"; output; 
macrovariable="cuminc"; value="&cuminc"; remark="Requests cumulative incidence curves"; output; remark="";
macrovariable="action"; value="&action"; 
	%if %upcase(%substr(&action,1,1))=C %then %do; remark="No output produced (see log file)"; %end; 
	%if %upcase(%substr(&action,1,1))=E %then %do; remark="Fine-Gray model computed."; %end; output; remark="";
macrovariable="weights"; value="&weights"; 
	if &weights=0 then remark="Standard model, no weighting of risk sets"; 
	if &weights=1 then remark="IPC weighting of risk sets"; 
	if &weights=2 then remark="AHR weighting of risk sets"; 
	output; remark="";
macrovariable="clean"; value="&clean"; 
	if &clean=1 then remark="Unnecessary variables removed"; 
	if &clean=0 then remark="Unnecessary variables not removed"; 
	output;
macrovariable="call"; value="&call"; remark="Data set with this call`s macro options"; output; remark="";
macrovariable="out"; value="&out"; remark="Output data set for standard Fine-Gray model"; output; remark="";
%if &weights ne 0 %then %do;
 macrovariable=""; value="&out._w"; remark="Output data set with weighted risk sets"; output; remark="";
%end;
macrovariable="missing"; value="&missing"; 
	%if %upcase(%substr(&missing,1,1))=K %then %do; remark="Keep observations with missing covariate values"; %end; 
	%if %upcase(%substr(&missing,1,1))=D %then %do; remark="Delete observations with missing covariate values"; %end;
	output; remark="";
macrovariable="statustab"; value="&statustab"; 
	if &statustab then remark="Summary of status variable requested"; 
	else remark="No summary of status variable requested"; 
	output; remark="";
macrovariable="delwork"; value="&delwork"; 
	if &statustab then remark="Temporary data sets deleted on exit"; 
	else remark="Temporary data sets available on exit"; 
	output; remark="";
macrovariable="--------"; 
value="------------"; 
remark="-------------------------------------------"; 
output; 
remark="";

macrovariable="macro version"; value="&version"; output;
macrovariable="      build"; value="&build"; output;
run;


proc print noobs label width=uniform;
title3 "The PSHREG macro: summary of macro options";
var macrovariable value remark;
run;
title3;





* retranslating the new option "action" into the old "finegray" which is used in the remainder of the code;

%if %upcase(%substr(&action,1,1)) = C %then %let finegray=0;
%if %upcase(%substr(&action,1,1)) = E %then %let finegray=1;


* setting firth to 1 will include the keyword FIRTH as option in the model statement of PROC PHREG;

%if &firth=1 %then %let options=%str(&options firth);


* if there were any options specified, add a slash and 
  then the string can be directly passed to PROC PHREGs model statement;

%if &options ne  %then %let options=%str(/ &options);


* the weights option is experimental, it allows to calculate time- or population-averaged 
  subdistribution hazard ratios;

%let ipcw=0;
%let ahr=0;
%if &weights=1 %then %let ipcw=1;
%if &weights=2 %then %let ahr=1;
%if &ahr=1 %then %let ipcw=1;


* now a copy of the data set is created and called _work. It contains all relevant information 
  needed by the macro.;

data _work;
set &data;
if _n_=1 then do; _cummissy_=0; _cummissx_=0; end;
_misslinex_=0;  ** this is to count missing values in the covariates;
_missliney_=0;  ** this is to count missing values in the outcome variable;
file log;
if &time = . then do;
	_missliney_=1;
end;
if &cens = . then do;
	_missliney_=1;
end;
if _cummissy_=0 and _missliney_=1 then put "WARNING: missing values in time or status variable.";
%if %upcase(%substr(&missing,1,1))=D %then %do;
	%do j=1 %to &nvar;
		if &&var&j=. then do;
			_misslinex_=1;
		end;
	%end;
	if _cummissx_=0 and _misslinex_=1 then put "WARNING: missing values in covariates.";
%end;

_cummissx_+_misslinex_;
_cummissy_+_missliney_;
_cummiss_+(1-(1-_misslinex_)*(1-_missliney_));
call symput("cummissx",_cummissx_);  * output the sum of missing values as macro variables;
call symput("cummissy",_cummissy_);
call symput("cummiss",_cummiss_);
* "clean" will only keep variables that were mentioned in the macro call;
%if &clean %then %do; 
	keep &time &cens &varlist &cengroup &by &id _misslinex_ _cummissx_ _missliney_ _cummissy_ _cummiss_ &admin; 
%end;
run;

%put NOTE: &cummissx lines deleted because of missing values in the independent variables.;
%put NOTE: &cummissy lines deleted because of missing values in the outcome variables.;

data _work;
set _work;
if _misslinex_ ne 1 and _missliney_ ne 1;      ** delete lines with missing values;
*drop _misslinex_ _missliney_ _cummissy_ _cummissx_ _cummiss_;
run;





data _work;
set _work;
if &cens=&failcode then _status=1;                   ** generate new status variable with 1 for a failure and 0 for censored observation;
else if &cens=&cencode then _status=0;
 else _status=2;                                     ** and 2 for competing event;
_time=&time;
%if &admin ne %then %do; if _status=2 then _time=&admin; %end;
%if &id = %then %do;
 %let ID=_id_;
 _id_=_n_;
%end;
%if &by = %then %do;
 %let byempty=1;                                    ** "remember" if by variable was specified;
 %let by=_by_;                                      ** if not, generate an artificial constant BY variable;
 _by_=1;
%end;
%else %do;
 %let byempty=0;
%end;
run;

%if &byempty ne 1 %then %do;
	proc sort data=_work out=_work;                 ** if a BY variable was specified, sort data by this variable 
														(thanks to Lovedeep Gondara for this remark);
	by &by;
	run;
	proc freq data=_work;
	tables &by*_status/noprint out=_censtab;        ** generate a table with frequencies of the status indicator;
	run;
%end;
%else %do;
	proc freq data=_work;
	tables _status/noprint out=_censtab;            ** generate a table with frequencies of the status indicator;
	run;
%end;

%if %substr(%upcase(&tiedcens),1,1)=A %then %do;    ** if tiedcens=after, then censored observations should be interpreted 
														as occuring after  tied event times;
	data _workt;                                    ** this is achieved by first computing the minimal non-zero difference in the time-variable,;
	set _work;                                      ** and then adding one tenth of that value to the censored times.;
	keep &by _time;
	run;
	proc sort data=_workt out=_workt;
	by &by _time;
	run;
	data _workt;
	set _workt;
	by &by;
	__diff=_time-lag(_time);
	if first.&by then __diff=0;
	run;
	proc means noprint;
	var __diff;
	output out=_mindiff min=__diff;
	where __diff ne 0;
	run;
	data _mindiff;
	set _mindiff;
	call symput("MINDIFF", __diff);
	run;
	data _work;
	set _work;
	if _status=0 then _time=_time+&mindiff/10;
	run;
%end;


%let warnadmin=0;
%if &admin ne %then %do;
	%if &ipcw = 1 | &ahr =1 %then %do;
		%let ipcw=0;
		%let ahr=0;
		%put WARNING: Weighting option ignored because of ADMIN option.;
	%end;
	data &out;
	set _work;
	if _n_=1 then warn=0;
	_weight_=1;
	_start_ =0;
	_stop_= _time;
	&censcrr=_status;
	if _status=2 then &censcrr=0;
	if _status=0 AND &time ne &admin then warn=1;
	retain warn;
	call symput("warnadmin",warn);
	run;

%end;
%else %do;

	**14-02-2014, JRA, add of this proc sort;
	proc sort data=_work;
	by &by &cengroup ;
	run;



	ods graphics off;
	proc lifetest data=_work noprint outsurv=_censdist;  ** compute the weights for the Fine-Gray model by Kaplan-Meier analysis with reversed;
	time _time*_status(1,2);                             ** meaning of the status indicator;
	by &by &cengroup;
	run;
	ods graphics;
	  

	data _censdist;	                                     ** generate a useable data set with the censoring distribution which can easily be;
	set _censdist;                                       ** merged with the original data;
	 by &by &cengroup _time;
	 if last._time;
	 if survival ne . then _Cdist_=survival;
	 retain _cdist_;
	 _time_=_time;
	keep _time_ _time _cdist_ &by &cengroup;
	 run;

	* this data set contains only events of interest and censored observations;
	data _dat01_crr;
	set _work;
	&censcrr=_status;
	if _status=0 or _status=1;
	_start_=0;
	_stop_=&time;
	*intervalcens=cens;
	run;

	* this data set contains the competing events;
	data _dat2_crr;
	set _work;
	if _status=2;
	run;
	proc sort;
	by &by _time;
	run;

	**14-02-2014, JRA, modif of this proc sort;
	proc sort;
	by &by &cengroup _time ;
	run;


	data _dat2_crr;
	merge _dat2_crr _censdist;
	by &by &cengroup _time;
	_denom_w_=_cdist_;
	if _status=2;
	drop _time_ _cdist_;
	run;
	proc sort;
	by &by &cengroup _time;
	run;

	data _censdist2;
	set _censdist;
	drop _time;
	run;

	* now we need proc sql to make an efficient merging of the data sets (assigning the weights to the original data);

	**14-02-2014, JRA, add of &cengroup in the select;

	PROC SQL;
	create table _censdist2_tmp as
	select &by, _cdist_,  %if &cengroup ne %then %do; &cengroup, %end;
	  min(_time_) as _time_min_,
	  max(_time_) as _time_max_
	from _censdist2
	group by &by, %if &cengroup ne %then %do; &cengroup, %end; _cdist_
	order by &by, %if &cengroup ne %then %do; &cengroup, %end; _time_max_ asc
	;

	data _censdist2_tmp;
	set _censdist2_tmp;
	_line_=_n_;
	ruN;

	proc sql;
	create table _censdist2_1_tmp as
	select t1.&by, %if &cengroup ne %then %do; t1.&cengroup, %end; t1._cdist_, t2._time_max_ as _time_min_, t1._time_max_, t1._line_
	from _censdist2_tmp t1
	left outer join _censdist2_tmp t2
	  on t1._line_ = t2._line_ + 1 and t1.&by = t2.&by %if &cengroup ne %then %do; and t1.&cengroup=t2.&cengroup %end;
	;

	create table _dat2a_crr as
	 SELECT d.*, c._cdist_, c._time_min_, c._time_max_ FROM _dat2_crr d, _censdist2_1_tmp c
	 where d._time <= c._time_max_ and d.&by = c.&by %if &cengroup ne %then %do; and d.&cengroup=c.&cengroup %end;
	;


	proc sort data=_dat2a_crr out=_dat2a_crr;
	by &by &id _time_max_;
	run;


	* the core of the macro: generating multiple lines with decreasing weights after each original competing event;
	data _dat2a_crr;
	set _dat2a_crr;
	_weight_=_cdist_/_denom_w_;
	&censcrr=0;
	if _time_min_<_time then do;
	 _start_=0;
	 _stop_=_time;  *** vorsicht, _start_ auf 0 setzen!***;
	 output;
	 _start_=_time;
	 _stop_=_time_max_;
	 output;
	end;
	else do;
	  _start_=_time_min_;
	  _stop_=_time_max_;
	  output;
	end;
	run;

	data &out;
	set _dat01_crr _dat2a_crr;
	if _weight_=. then _weight_=1;
	if _start_ < _stop_;
	drop _time_min_ _time_max_ _denom_w_ _cdist_ _status _time;
	run;

	proc sort;
	by &by;
	run;
%end;

%if &cummiss ne 0 %then %do;   * Report on missing values;
	data __x;
	cummissx=&cummissx;
	cummissy=&cummissy;
	cummiss=&cummiss;
	remark="Observations deleted in input data set because of missing values:";
	label cummissx="Covariates" cummissy="Outcome (&time, &cens)" cummiss="Total";
	run;
	proc print data=__x noobs label;
	%if %upcase(%substr(&missing,1,1))=D %then %do;
		title3 "The PSHREG macro: Summary of missing values in covariates and outcome variables";
		var remark cummissx cummissy cummiss;
	%end;
	%else %do;
		title3 "The PSHREG macro: Summary of missing values in outcome variables";
		var remark cummiss;
	%end;
	run;
%end;

%if &statustab=1 %then %do;  * report on status indicator;
	proc format;
	value status 0="Censored" 1="Events of interest" 2="Competing events";
	run;
	proc print data=_censtab;
	var %if &byempty ne 1 %then %do; &by %end; _status count percent;
	format _status status.;
	title3 "The PSHREG macro: Summary of status variable";
	label _status="&cens";
	run;
%end;

%let x=	%scan(&varlist,1);    * x is the first variable in the variable list. 
                                it will be used for nonparametric cumulative incidence curve, if needed.;

%if &cuminc=1 %then %do;     * generate nonparametric eventtype-specific cumulative incidence curves;
		proc phreg data=&out noprint;
		model (_start_,_stop_)*&censcrr(0)=;
		weight _weight_;
		strata &x;
		%if &byempty=0 %then %do; by &by; %end;
		baseline out=_cuminc survival=_surv /method=EMP;
		run;

		data _cuminc;
		set _cuminc;
		_cuminc_=1-_surv;
		&time=_stop_;
		label _cuminc_="Cumulative incidence";
		drop _stop_ _surv;
		run;

		symbol1 i=steplj line=1 c=black;
		symbol2 i=steplj line=2 c=black;
		proc gplot data=_cuminc;
		plot _cuminc_*&time=&x;        * plot the cumulative incidence curves;
		%if &byempty=0 %then %do; by &by; %end;
		title3 "The PSHREG macro: Cumulative incidence curves";
		run;

%end;


%if &ipcw %then %do;    * this is for the experimental time-averaged analysis;
		proc sql;
		*** selects all event-of-interest times ***;
		create table _events as
		select &by, %if &cengroup ne %then %do; &cengroup, %end; _time as _time, count(*) as _cnt from _work
		where _status=1
		group by &by, %if &cengroup ne %then %do; &cengroup, %end; _time
		order by &by, %if &cengroup ne %then %do; &cengroup, %end; _time
		;

		*** adds ipcw _weight_s to event times ***;
		create table _ipcw as
		 select t.*, 1./c._cdist_ as ipcw
		 from _events t, _censdist c
  	     where t.&by = c.&by %if &cengroup ne %then %do; and t.&cengroup=c.&cengroup %end; and t._time = c._time_
		 ;

		*** blows up data set, stratifies by risk set*;
		create table &out._w as
		 select d.*, e._time as _riskset_, e.ipcw*d._weight_ as _ipcweight_,(d._stop_=e._time and d.&censcrr=1) as _wcens_    
         from &out d
		 inner join _ipcw e
		 on d.&by = e.&by %if &cengroup ne %then %do; and d.&cengroup = e.&cengroup %end; and d._start_ < e._time and d._stop_>=e._time
		 order by d.&by, %if &cengroup ne %then %do; d.&cengroup, %end; e._time, d._stop_
		 ;
		

%end;


%if &ahr=1 %then %do;    * this is for the experimental population-averaged analysis. Shh, do you want to know a secret? 
						  You can use that option to fit a weighted Cox model (Schemper, Wakounig, Heinze, StatMed 2009)
                          in case of no competing events;
		proc phreg data=&out noprint;
		model (_start_,_stop_)*&censcrr(0)=;
		weight _weight_;
		baseline out=_smarg survival=_surv /method=EMP;
		by &by; 
		run;
		
		data _smarg;
		set _smarg;
		_line_=_n_;
		ruN;

		proc sql;
		* creates left-continuous pseudo-survival function estimate;
		create table _smarg_tmp as
		select t1.&by, t1._surv, t2._stop_ as _stop_, t1._stop_ as _stop_old, t1._line_
		from _smarg t1
		inner join _smarg t2
		  on t1._line_ = t2._line_ - 1 and t1.&by = t2.&by
		;



		data &out._w_tmp;
		set &out._w;
		run;

		proc sql;
		create table &out._w as
		select d.*, e._surv as _sweight_, e._surv*d._ipcweight_ as _ahrweight_ from &out._w_tmp d
        join _smarg_tmp e
		on d.&by=e.&by and e._stop_=d._riskset_
		order by d.&by, d._riskset_
		;
		drop table &out._w_tmp
		;
		** some notes from the program developer in his own language:;
        ***gewichte aus left-continuous Schätzer nehmen und mit den Gewichten für G-1 vermengen***;
		*** datensatz aufblasen sodass es für jedenevent of interest ein eigenes Stratum gibt wo alle Pat at risk drinnen sind;
%end;







%if &finegray=1 %then %do;   ** now the model is actually calculated;
	%if &ipcw=0 %then %do;
		 proc phreg data=&out %if &firth=0 %then %do; covs(aggregate) %end;;
		 title3 "The PSHREG macro: Fine-Gray model";
		 %if &class ne %then %do; class &class; %end;
		 model (_start_,_stop_)*&censcrr(0)=&varlist &options;
		 weight _weight_;
		 id &id;
		 %if &byempty ne 1 %then %do; by &by; %end;
 		 run;
	%end;
	%if &ipcw=1 & &ahr=0 %then %do;   *** time-averaged analysis;
		 proc phreg data=&out._w %if &firth=0 %then %do; covs(aggregate) %end;;
		 title3 "The PSHREG macro: Fine-Gray model with IPC weighting of risksets";
		 %if &class ne %then %do; class &class; %end;
		 model (_start_,_stop_)*_wcens_(0)=&varlist &options;
		 weight _ipcweight_;
		 strata &cengroup _riskset_;
		 id &id;
 		 %if &byempty ne 1 %then %do; by &by; %end;
		 run;
	%end;
	%if &ahr=1 %then %do;            *** population-averaged analysis;
		 proc phreg data=&out._w %if &firth=0 %then %do; covs(aggregate) %end;;
		 title3 "The PSHREG macro: Fine-Gray model with AHR weighting of risksets";
		 %if &class ne %then %do; class &class; %end;
		 model (_start_,_stop_)*_wcens_(0)=&varlist &options;
		 weight _ahrweight_;
		 strata &cengroup _riskset_;
		 id &id;
		 %if &byempty ne 1 %then %do; by &by; %end;
		 run;
	%end;
%end;

%if &delwork=1 %then %do;
	* garbage collection;
	proc datasets nolist; ***thanks to Akos F. Pap, Wuppertal, D, for indicating that nolist is the proper option;
	delete _dat01_crr _dat2a_crr _censdist _censdist2 _censdist2_1_tmp _censdist2_tmp _censtab _dat2_crr _events _work __x ;
	%if &weights ne 0 %then %do;
	 delete _ipcw;
	%end;
	%if &weights = 2 %then %do;
	 delete _smarg _smarg_tmp;
	%end;
	run;
	quit;
%end;


*** if cuminc=0 then we tell the user (in the Log window) how to obtain the nonparametric cumulative incidence curves;
%if &cuminc=0 %then %do;
	data __x;
	file log;
	put "NOTE: Code for cumulative incidence curves stratified by &x :";
	put "NOTE: (This code will only work if &x is a categorical variable.)";
	put "PROC PHREG DATA=&out ;";
		put "MODEL (_start_,_stop_)*&censcrr(0)=;";
		put "WEIGHT _weight_;";
		put "STRATA &x; *** check if &x is categorical! ***;";
		put "BASELINE OUT=_cuminc SURVIVAL=_surv /METHOD=EMP;";
		%if &byempty=0 %then %do; put "BY &by;"; %end;
		put "RUN;";

		put "DATA _cuminc;";
		put "SET _cuminc;";
		put "_cuminc_=1-_surv;";
		put "&time=_stop_;";
		put "LABEL _cuminc_=""Cumulative incidence"";";
		put "DROP _stop_ _surv;";
		put "RUN;";

		put "SYMBOL1 I=steplj LINE=1 C=black;";
		put "SYMBOL2 I=steplj LINE=2 C=black;";
		put "PROC GPLOT DATA=_cuminc;";
		put "PLOT _cuminc_*&time=&x;";
		%if &byempty=0 %then %do; put "BY &by;"; %end;
		put "RUN;";
	run;


%end;

%if &finegray=0 %then  %do;     *** if action=code then we output the code into the Log window;
	data __x;
	file log;
	put "NOTE: for the Fine-Gray model: call phreg by:";
	%if &firth=0 %then %do; put "PROC PHREG DATA=&out COVS(AGGREGATE);"; %end;
	%else %do; put "PROC PHREG DATA=&out ;"; %end;
	%if &class ne %then %do; put "CLASS &class;"; %end;
	put "MODEL (_start_,_stop_)*&censcrr(0)=&varlist &options;";
	put "ID &id;";
	put "WEIGHT _weight_;";
	%if &byempty ne 1 %then %do; put "BY &by;"; %end;
	put "RUN;";
	%if &ipcw=1 %then %do;
		put "NOTE: for the weighted Fine-Gray model with inverse probability of censoring weights: call phreg by:";
		%if &firth=0 %then %do; put "PROC PHREG DATA=&out._w COVS(AGGREGATE);"; %end;
		%else %do; put "PROC PHREG DATA=&out._w ;"; %end;
 	    %if &class ne %then %do; put "CLASS &class;"; %end;
		put "MODEL (_start_,_stop_)*_wcens_(0)=&varlist &options;";
		put "STRATA _riskset_;";
		put "ID &id;";
		put "WEIGHT _ipcweight_;";
		%if &byempty ne 1 %then %do; put "BY &by;"; %end;
		put "RUN;";
		%if &ahr=1 %then %do;
			put "NOTE: for the weighted Fine-Gray model to compute average subdistribution hazard ratios: call phreg by:";
			%if &firth=0 %then %do; put "PROC PHREG DATA=&out._w COVS(AGGREGATE);"; %end;
			%else %do; put "PROC PHREG DATA=&out._w ;"; %end;
			%if &class ne %then %do; put "CLASS &class;"; %end;
            put "MODEL (_start_,_stop_)*_wcens_(0)=&varlist &options;";
			put "STRATA _riskset_;";
			put "ID &id;";
			put "WEIGHT _ahrweight_;";
			%if &byempty ne 1 %then %do; put "BY &by;"; %end;
			put "RUN;";
		%end;
	%end;

	run;
%end;

%if &warnadmin = 1 %then %do;
		%put WARNING: some censored times found which are not equal to administrative censoring times. Random censoring? ;
%end;


proc datasets nolist;   ***thanks to Akos F. Pap, Wuppertal, D, for indicating that nolist is the proper option;
delete __x; ** die letzte Kuh sperrts Türl zu...;  *** delete the last pseudo data set;
run;
quit;
title3;
%mend;

