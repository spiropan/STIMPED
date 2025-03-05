%% Download datasets from the below links
% COVID-19_vaccinations_sorted downloaded from
% https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc

% Excess deaths downloaded from 
% https://data.cdc.gov/NCHS/AH-Excess-Deaths-by-Sex-Age-and-Race-and-Hispanic-/m74n-4hbs/about_data

% Daily cases per 1M population downloaded from
% https://ourworldindata.org/covid-cases

%% Grab excess mortality (deaths) for all ages across 2020-2022
clear

% This file from https://data.cdc.gov/NCHS/AH-Excess-Deaths-by-Sex-Age-and-Race-and-Hispanic-/m74n-4hbs/about_data
file='./Data/AH_Excess_Deaths_by_Sex__Age__and_Race_and_Hispanic_Origin_20250131.csv';
dat=readtable(file);

% Define years
years=[2020,2021,2022];

% Grab age groups
ages=unique(dat.AgeGroup); ages(strcmp(ages,'Not stated'),:)=[];

% Initialize 3D array that will hold excess deaths 
% First dimenions rows=week, 2nd dim columns=age group
% 3rd dimension = both sexes, males and females
EM=[];

% Initialize array that will hold the week labels
EM_week=[];

% Loop over all years
for y=1:length(years)
    disp(['Grabbing excess deaths for year: ' int2str(years(y))])
    
    % Grab the subset of data in a given year
    subdat=dat(dat.MMWRyear == years(y),:); 
    weeks=unique(subdat.MMWRweek);
    
    % Make a week X age tmp array of excess deaths
    tmp_EM=zeros(size(weeks,1),size(ages,1),3);
    tmp_week=unique(subdat.Weekending);
    
    % Loop over all weeks in Year 2021
    for w=1:length(weeks) % exclude ages 'not stated'
        
        % Now loop over all ages    
        for a=1:length(ages)
            tmp_EM(w,a,1)=table2array(subdat(subdat.MMWRweek==weeks(w) & strcmp(subdat.AgeGroup,ages{a})...
                & strcmp(subdat.RaceEthnicity,'All Race/Ethnicity Groups') & strcmp(subdat.Sex,'All Sexes'),...
                {'NumberAboveAverage_weighted_'})); % {'PercentAboveAverage_unweighted_'}));
            tmp_EM(w,a,2)=table2array(subdat(subdat.MMWRweek==weeks(w) & strcmp(subdat.AgeGroup,ages{a})...
                & strcmp(subdat.RaceEthnicity,'All Race/Ethnicity Groups') & strcmp(subdat.Sex,'Male (M)'),...
                {'NumberAboveAverage_weighted_'})); % {'PercentAboveAverage_unweighted_'}));
            tmp_EM(w,a,3)=table2array(subdat(subdat.MMWRweek==weeks(w) & strcmp(subdat.AgeGroup,ages{a})...
                & strcmp(subdat.RaceEthnicity,'All Race/Ethnicity Groups') & strcmp(subdat.Sex,'Female (F)'),...
                {'NumberAboveAverage_weighted_'})); % {'PercentAboveAverage_unweighted_'}));
        end
    end
    
    % Concatenate arrays across each year
    EM=[EM; tmp_EM];
    EM_week=[EM_week; tmp_week]; 
end

%% Here grab weekly vaccinations for 2021-2022
vax=readtable('./Data/COVID-19_Vaccinations_in_the_United_States_Jurisdiction_20250123.csv');

% Grab the vaccination dates and sort them
vax_dates=unique(vax.Date);
[Y,I]=sort(datenum(vax_dates));
vax_dates=vax_dates(I);

% CDC stopped recording daily vaccinations from 6/16/2022 onwards. This
% cell records the dates on which they recorded the vaccinations so that
% they can be added manually to the vaccination weeks vector
SPECIAL={'6/15/2022';'6/22/2022';'6/29/2022';'7/6/2022';'7/13/2022';...
            '7/20/2022';'7/27/2022';'8/3/2022';'8/10/2022';'8/17/2022';'8/24/2022';...
            '8/31/2022';'9/7/2022';'9/14/2022';'9/21/2022';'9/28/2022';'10/5/2022';...
            '10/12/2022';'10/19/2022';'10/26/2022';'11/2/2022';'11/9/2022';'11/16/2022';...
            '11/23/2022';'11/30/2022';'12/7/2022';'12/14/2022';'12/21/2022';'12/28/2022'};
        
vax_weeks=[EM_week(1:129,1); SPECIAL(1:end)]; vax_weeks=intersect(vax_dates,vax_weeks);

vax_counts=[];
% For each week in vax_weeks, convert the cumulative Administered doses to
% the total given in each week

for w=1:length(vax_weeks)
   
    % ----- Administered (total) ----- %
    if w==1
        cur_week = table2array(vax(strcmp(vax.Location,'US') & (vax.Date == vax_weeks(w)),...
        {'Administered','Administered_5Plus','Administered_12Plus','Administered_18Plus','Administered_65Plus'}));
         cur_week(isnan(cur_week))=0;
         vax_counts(w,:)=cur_week;
    else
        cur_week = table2array(vax(strcmp(vax.Location,'US') & (vax.Date == vax_weeks(w)),...
            {'Administered','Administered_5Plus','Administered_12Plus','Administered_18Plus','Administered_65Plus'}));
        cur_week(isnan(cur_week))=0;
        
        pre_week = table2array(vax(strcmp(vax.Location,'US') & (vax.Date == vax_weeks(w-1)),...
            {'Administered','Administered_5Plus','Administered_12Plus','Administered_18Plus','Administered_65Plus'}));
        pre_week(isnan(pre_week))=0;
        vax_counts(w,:) = cur_week-pre_week;
    end
end

% Zero pad vax_counts to match the EM_weeks
padsize=size(EM,1)-size(vax_counts,1); vax_counts=padarray(vax_counts,padsize,0,'pre');

%% Here grab weekly COVID cases for 2020-2022
cas=readtable('./Data/daily-new-confirmed-covid-19-cases-per-million-people.csv');

% Can alternativelty use the weekly (non 7-day rolling average) cases
%cas=readtable('./Data/weekly-confirmed-covid-19-cases.csv');

% Grab the dates from the spreadsheet
cas_dates=unique(cas.Day);

% Apply multiple lags from 14 days to no lag and make cas_counts a 157x14 matrix where each column
% stores the values for each lag
lags = [14:-1:0];
cas_counts=lagmatrix(cas.DailyNewConfirmedCasesOfCOVID_19PerMillionPeople_rolling7_dayAv,lags);

% Determine which days for the cases correspond to weeks that are included
% in the excess deaths vector
cas_weeks=intersect(cas_dates,EM_week);
cas_idx=ismember(cas.Day,cas_weeks);

% Below if using weekly cases instead 
%cas_counts=cas.WeeklyCases(cas_idx);

% Now downsample the lagged COVID cases to a weekly resolution
cas_counts=cas_counts(cas_idx,:);

% Remove any NaN
cas_counts(isnan(cas_counts))=0;

% zero pad the arrays to make the same size as EM
padsize=size(EM,1)-size(cas_counts,1); cas_counts=padarray(cas_counts,padsize,0,'pre');

%% Plot all variables (zscored)
% Plot excess deaths for all ages (all ages is the indexed by 17 in our
% array). 
plot(smooth(zscore(EM(:,17,1)))); hold on; plot(smooth(zscore(vax_counts(:,1))));
hold on; plot(smooth(zscore(cas_counts(:,1))));

%% Here determine the best fitting day lag for cases with single COVID case regressor
% Fit for each case lag and store the variance explained

% Array to store R^2 for each model
ve=[]; 

% Struct to store the model outputs from fitlm for each lagged version of
% the covid cases vector
dat=struct();

% Loop over all lags, fit a model, and store the results
for i=1:length(lags)
    DT = [table(cas_counts(:,i),EM(:,17,1),'VariableNames',{'Cases','EM'})]; % Table
    M0=fitlm(DT,'Intercept',true); 
    dat(i).M0=M0;
    ve(i,1)=M0.Rsquared.Ordinary;
end

% Optimal lag in days. This yields 4 days.
optimal_lag=lags(find(ve==max(ve)));

% The optimal lag index will be used in all subsequent analyses
ol_idx=find(ve==max(ve));

%% Display results from the first model (M0) that includes just a single COVID regressor
dat(ol_idx).M0

%% Improve the model by including a separate terms of each COVID wave
% The below points were determined manually by selecting the local minima
% based on a plot of the (non-lagged) cas_counts
% plot(cas_counts(:,15));
inflections=[1,24,37,63,79,97,119,148,158];

cas_waves=zeros(size(cas_counts,1),8,size(cas_counts,2)); % make 3D array, with lag as 3rd dimensions
for ii=1:length(lags)
    for i=1:length(inflections)-1
        cas_waves(inflections(i):inflections(i+1)-1,i,ii)=cas_counts(inflections(i):inflections(i+1)-1,ii);
    end
end

% Set the age index to "all ages" and optimal lag index to i
a=17; i=ol_idx;

% Create data table for the M1 with 8 seperate COVID waves, with a 4 day
% lag
cas_labels={'w1','w2','w3','w4','w5','w6','w7','w8','EM'};
DT1 = [table(cas_waves(:,1,i),cas_waves(:,2,i),cas_waves(:,3,i),cas_waves(:,4,i),...
       cas_waves(:,5,i),cas_waves(:,6,i),cas_waves(:,7,i),cas_waves(:,8,i),EM(:,17,1),...
       'VariableNames',cas_labels)]; 
   
M1=fitlm(DT1)

%% Now estimate the model with the vaccine term added   
dos_labels={'w1','w2','w3','w4','w5','w6','w7','w8','doses','EM'};
DT2 = [table(cas_waves(:,1,i),cas_waves(:,2,i),cas_waves(:,3,i),cas_waves(:,4,i),...
        cas_waves(:,5,i),cas_waves(:,6,i),cas_waves(:,7,i),cas_waves(:,8,i),vax_counts(:,1),EM(:,17,1),...
        'VariableNames',dos_labels)]; % Table

M2=fitlm(DT2)

%% Estimate excess deaths from full model, all ages and both sexes (M1 above)
% Estimate "COVID deaths" in model with a single COVID regressor
M0=dat(ol_idx).M0;
beta=table2array(M0.Coefficients(2,1));
cas_vec=table2array(M0.Variables(:,1)); cas_vec(isnan(cas_vec))=[];
total_cases=sum(cas_vec)*7*331;
covid_deaths=beta*total_cases/1000000*1036 % 156 is days through June 6, 2020
cfr=covid_deaths/total_cases;

% Estimate excess deaths within a week post-injection
vfr=M2.Coefficients.Estimate(end)
doses=table2array(M2.Variables(:,9));
deaths=vfr*sum(doses)

%% Plot the covid waves and their labels to aid in interpretation
dt = datetime(string(EM_week),'Format','dd-MM-yy'); 
fig=figure;
ax=axes(fig);
X=cas_waves(:,:,15);
varNames={'w1','w2','w3','w4','w5','w6','w7','w8'};
plot(dt,[X],'LineWidth',2);
set(ax,'XTick',EM_week(1:16:end));
datetick(ax,'x','mm-yyyy','keepticks');
hold('on');
xlabel('Date')
ylabel('COVID cases per 1M population per day')
legend([varNames],'Location','NW')
title(['\bf COVID wave regressors'])
axis('tight')
grid('on')

%% Plot just the outcome variable
dt = datetime(string(EM_week),'Format','dd-MM-yy'); 
fig=figure;
ax=axes(fig);
X=EM(:,17,1);
plot(dt,[X],'LineWidth',2);
set(ax,'XTick',EM_week(1:16:end));
datetick(ax,'x','mm-yyyy','keepticks');
hold('on');
xlabel('Date')
ylabel('Weekly excess deaths')
title(['\bf Weekly excess deaths (vs. 2015-2019 avgs, both sexes, all ages)'])
axis('tight')
grid('on')

%% Plot excess deaths, COVID waves and vaccine doses all on the same plot
dt = datetime(string(EM_week),'Format','dd-MM-yy'); 
fig=figure;
ax=axes(fig);
X=zscore([EM(:,17,1) cas_counts(:,15) vax_counts(:,1)]);
Names={'Weekly excess deaths: both sexes, all ages','COVID cases per 1M per day','Doses administered'}
plot(dt,[X],'LineWidth',2);
set(ax,'XTick',EM_week(1:16:end));
datetick(ax,'x','mm-yyyy','keepticks');
hold('on');
xlabel('Date')
ylabel('Zscores')
legend([varNames],'Location','NW')
title(['\bf Weekly excess deaths and COVID cases per 1M per day'])
axis('tight')
grid('on')

%% Plot just w4
fig=figure;
ax=axes(fig);
starti=63; endi=79; % For Feb-June 2021
dt = datetime(string(EM_week(starti:endi)),'Format','dd-MM-yy');

% Zscored variables
Yz=zscore(EM(starti:endi,17,1));
Xz=zscore([cas_counts(starti:endi,15), vax_counts(starti:endi,1)])

% Raw variables
Y=EM(starti:endi,17,1);
X=[cas_counts(starti:endi,15), vax_counts(starti:endi,1)]
varNames={'Cases','Vax','EM'}
plot(dt,[Xz,Yz],'LineWidth',2);
hold('on');
xlabel('Date');
ylabel('Predictors and excess mortality (zscores)');
legend([varNames],'Location','NW')
grid('on')

% Interrogate the time series chunks further
% Look at correlation between Vax and Cases
[r,p]=corr(X(:,1),X(:,2),'type','spearman')

% Estimate fit just within w4
[b,dev,stats]=glmfit([X],[Y])

% Derive estimates of VFR just using difference in rise in excess deaths
% divided by rise in vaccine doses during w4
(max(smooth(Y))-min(smooth(Y)))/(max(smooth(X(:,2))-min(smooth(X(1,2)))))

(max(Y)-min(Y))/(max(X(:,2)-X(1,2)))

%% Create table of vaccine coefficients by age and lag
% Lag vaccine term from 0-10 weeks
vlags=[0:10]; 
% index of the optimum lag for COVID cases (4 days)
i=ol_idx; 
clear clms hacs

% In addition to OLS, we will estimate heteroscedasticity and 
% autocorrelation consistent Newey-West covariance matrices as outlined 
% here: https://www.mathworks.com/help/econ/hac.html#btrsadq-2

% Compute maxlag bandwith as in [2] Andrews, D. W. K., and J. C. Monohan.
% "An Improved Heteroskedasticity and Autocorrelation Consistent 
% Covariance Matrix Estimator." Econometrica. Vol. 60, 1992, pp. 953–966.
% as cited in https://www.mathworks.com/help/econ/hac.html#btrsadq-2
T=height(DT);
maxLag = floor(4*(T/100)^(2/9));

% initialize array of vaccine coefficients and standard errors
clms.coeff=zeros(length(ages),length(vlags)); 
clms.se=zeros(length(ages),length(vlags));
hacs.coeff=zeros(length(ages),length(vlags)); 
hacs.se=zeros(length(ages),length(vlags));

% initialize column of table labels for age
for v=1:length(vlags)
    for a=1:length(ages)
        disp(['Working on lag: ' int2str(v) ', age ' int2str(a)])
        vax_lag=lagmatrix(vax_counts(:,1),vlags(v));
        dos_labels={'w1','w2','w3','w4','w5','w6','w7','w8','doses','EM'};
        DT = [table(cas_waves(:,1,i),cas_waves(:,2,i),cas_waves(:,3,i),cas_waves(:,4,i),...
            cas_waves(:,5,i),cas_waves(:,6,i),cas_waves(:,7,i),cas_waves(:,8,i),vax_lag,EM(:,a,1),...
            'VariableNames',dos_labels)]; % Table
        M=fitlm(DT);
        [EstCoeffCov,se,coeff]=hac(DT,'Bandwidth',maxLag+1);
        hacs.se(a,v)=se(end); hacs.coeff(a,v)=coeff(end);
        clms.coeff(a,v)=M.Coefficients.Estimate(end);
        clms.se(a,v)=M.Coefficients.SE(end);
        clm_tmp=[num2str(clms.coeff(a,v),'%0.1d') ' (' num2str(clms.se(a,v),'%0.1d') ')'];
        hac_tmp=[num2str(hacs.coeff(a,v),'%0.1d') ' (' num2str(hacs.se(a,v),'%0.1d') ')'];
        clms.tab{a,v}=clm_tmp;
        hacs.tab{a,v}=hac_tmp;
    end
end

for v=1:length(vlags)
    for a=1:length(ages)
        clms.tval(a,v)=clms.coeff(a,v)/clms.se(a,v);
        clms.pval(a,v)=2*(1-tcdf(abs(clms.tval(a,v)),147));
        clm_tmp=[num2str(clms.tval(a,v),'%0.2f') ' (' num2str(clms.pval(a,v),'%0.4f') ')'];
        clms2.tab{a,v}=clm_tmp;
        
        hacs.tval(a,v)=hacs.coeff(a,v)/hacs.se(a,v);
        hacs.pval(a,v)=2*(1-tcdf(abs(hacs.tval(a,v)),147));
        hac_tmp=[num2str(hacs.tval(a,v),'%0.2f') ' (' num2str(hacs.pval(a,v),'%0.4f') ')'];
        hacs2.tab{a,v}=hac_tmp;
    end
end
       
clms_tab=cell2table([ages clms2.tab],'VariableNames',{'Age','L0: tval (pval)','L1: tval (pval)','L2: tval (pval)',...
    'L3: tval (pval)','L4: tval (pval)','L5: tval (pval)','L6: tval (pval)','L7: tval (pval)','L8: tval (pval)',...
    'L9: tval (pval)','L10: tval (pval)'})

hacs_tab=cell2table([ages hacs2.tab],'VariableNames',{'Age','L0: tval (pval)','L1: tval (pval)','L2: tval (pval)',...
    'L3: tval (pval)','L4: tval (pval)','L5: tval (pval)','L6: tval (pval)','L7: tval (pval)','L8: tval (pval)',...
    'L9: tval (pval)','L10: tval (pval)'})

writetable(clms_tab,'CLM_Table.csv')
writetable(hacs_tab,'HAC_BW5_Table.csv')

[h crit_p hacs.adj_p]=fdr_bh(hacs.pval,0.05,'pdep');
[h crit_p clms.adj_p]=fdr_bh(clms.pval,0.05,'pdep');

% Write out HAC table with terms suriving p<0.05 FDR corrected 
for v=1:length(vlags)
    for a=1:length(ages)
        if hacs.adj_p(a,v) > 0.05
            hacs_tab(a,v+1)={'ns'};
        end
        if clms.adj_p(a,v) > 0.05
            clms_tab(a,v+1)={'ns'};
        end
    end
end

writetable(hacs_tab,'HAC_BW5_ADJ_Table.csv')
writetable(clms_tab,'CLM_ADJ_Table.csv')   

%% Plot residuals 
close_figures
DT = [table(cas_waves(:,1,i),cas_waves(:,2,i),cas_waves(:,3,i),cas_waves(:,4,i),...
    cas_waves(:,5,i),cas_waves(:,6,i),cas_waves(:,7,i),cas_waves(:,8,i),vax_counts(:,1),EM(:,17,1),...
    'VariableNames',dos_labels)]; % Table
M=fitlm(DT);

plotResiduals(M,"fitted")

tiledlayout(2,1)
nexttile
plotResiduals(M,"caseorder")
nexttile
plotResiduals(M,"lagged")

%% Write out data tables
dos_labels={'w1','w2','w3','w4','w5','w6','w7','w8','doses','EM'};
mkdir('Tables')
cd('./Tables')
for a=1:length(ages)
    
    DT = [table(cas_waves(:,1,i),cas_waves(:,2,i),cas_waves(:,3,i),cas_waves(:,4,i),...
    cas_waves(:,5,i),cas_waves(:,6,i),cas_waves(:,7,i),cas_waves(:,8,i),vax_counts(:,1),EM(:,a,1),...
    'VariableNames',dos_labels)]; % Table
    writetable(DT,['ages_' strrep(ages{a},' ','-') '-Data.csv']) 
end
cd('../')



        