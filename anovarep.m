function anovarep(x,varargin)
% ANOVAREP - Analysis of variance for repeated measures. 
% This function executes the analysis of variance when subjects underwent
% several treatments. This function is similar to ANOVA2 Matlab Function,
% but there are three differences:
% 1) the output of ANOVA table;
% 2) the graphical plot of anova table;
% 3) if p-value<alpha this function executes the Holm-Sidak test for multiple
% comparison test to highlight differences between treatments.
% 
% Syntax: 	ANOVAREP(X,ALPHA)
%      
%     Inputs:
%           X - data matrix. 
%           ALPHA - significance level (default = 0.05).
%     Outputs:
%           - Anova table.
%           - Graphical plot of anova table.
%           - Holm-Sidak test (eventually)
%
%      Example: anovarepdemo
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Anovarep: compute the Anova for repeated measures and
% Holm-Sidak test for multiple comparisons if Anova is positive. 
% http://www.mathworks.com/matlabcentral/fileexchange/18746

%Input error handling
p = inputParser;
addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d','real','finite','nonnan','nonempty'}));
addOptional(p,'alpha',0.05, @(x) validateattributes(x,{'numeric'},{'scalar','real','finite','nonnan','>',0,'<',1}));
parse(p,x,varargin{:});
alpha=p.Results.alpha;
clear p

[n,m]=size(x); %n=subjects; m=treatments;
T=mean(x); %Mean of treatment
S=mean(x,2); %Mean of subject
Xbar=sum(x(:))/(m*n); %General mean

SStra=m*sum((S-Xbar).^2); %Variability among subjects
GLtra=n-1; %degrees of freedom
MStra=SStra/GLtra; % variance of population estimed from subjects

SSentro=sum(sum((x-repmat(S,1,m)).^2,2)); %Variability within subjects
GLentro=n*(m-1); %degrees of freedom

SStot=SStra+SSentro; %Total Variability
GLtot=GLtra+GLentro; %Degrees of freedom

%The variability within subjects can be decomposed in:

SSt=n*sum((T-Xbar).^2); %Variability among treatments
GLt=m-1; %degrees of freedom
MSt=SSt/GLt; % variance of population estimed from treatments

SSr=SSentro-SSt; %Residual variability
GLr=(n-1)*(m-1); %degrees of freedom
MSr=SSr/GLr; % variance of population estimed from residuals

F=[MStra MSt]./MSr;
panova=1-fcdf(F,[GLtra GLt],GLr); %p-value

%Formatting for ANOVA Table printout.
rn={'Total';'Among subjects';'Within subjects';'Among treatments';'Residual'};
SS=[SStot; SStra; SSentro; SSt; SSr];
GL=[GLtot; GLtra; GLentro; GLt; GLr];
tanan=NaN(5,1);
MS=tanan; MS([2 4 5])=[MStra;MSt;MSr];
FA=tanan; FA([2 5])=[F(1);F(2)];
PV=tanan; PV([2 5])=[panova(1);panova(2)];
clear tanan
disp(table(rn,SS,GL,MS,FA,PV,'VariableNames',{'Variability','SS','df','MS','F','p_value'}))
if panova(2)<alpha
    disp('Almost one of the treatments is the cause of variation');
else
    disp('Variation is due to chance');
end

%Display bar a bar graph of anova table
y=[0 0 1 1];
hold on
h1=fill([0 SStot SStot 0],y,'c');
y=y+1;
h2=fill([0 SStra SStra 0],y,'y');
h3=fill([0 SSentro SSentro 0]+SStra,y,'r');
y=y+1;
h4=fill([0 SSt SSt 0]+SStra,y,'g');
h5=fill([0 SSr SSr 0]+SStra+SSt,y,'b');
hold off
legend([h1 h2 h3 h4 h5],'Total','Among subjects','Within subjects','Among treatments','Residual','Location','NorthWest')
title('Sources of Variability')
clear S Xbar SStra GLtra SSentro GLentro SStot GLtot SSt GLt MSt MStra SSr F rn SS GL MS FA PV
clear y h1 h2 h3 h4 h5 

%If Panova<alpha execute the Holm-Sidak test to close off the differences
%in anova
if panova(2)<alpha
    disp(' ')
    a=m-1; %rows of probability matrix
    c=0.5*m*(m-1); %max number of comparisons
    count=0; %counter
    p=ones(1,c); %preallocation of p-value vector
    pb{c,4} = []; %preallocation of p-value matrix
    denom=realsqrt(2*MSr/n); %the denominator is the same for each comparison
    for I=1:a
        for J=I+1:m
            count=count+1;
            t=abs(diff(T([I J])))/denom; %t-value
            p(count)=(1-tcdf(t,GLr))*2; %2-tailed p-value vector
            pb(count,1:2)={strcat(int2str(I),'-',int2str(J));p(count)}; %Matrix of the p-values
        end
    end
    clear a m t count I J %clear unnecessary variables
    [p,I]=sort(p); %sorted p-values
    pb=pb(I,:);
    J=1:c; %How many corrected alpha value?
    alphacorr=1-((1-alpha).^(1./(c-J+1))); %Sidak alpha corrected values
    %Compare the p-values with alpha corrected values. 
    %If p<a reject Ho; else don't reject Ho: no more comparisons are required.
    comp=1; %compare checker
    for J=1:c
        if comp %Comparison is required
            if p(J)<alphacorr(J)
                pb(J,3:4)={alphacorr(J);'Reject H0'};
            else
                pb(J,3:4)={alphacorr(J);'Fail to reject H0'};
                comp=0; %no more comparison are required
            end
        else %comparison is unnecessary
            pb(J,3:4)={'No comparison made';'H0 is accepted'};
        end
    end
    disp(cell2table(pb,'VariableNames',{'Comparison','p_value','Sidak_alpha','Comment'}))
end