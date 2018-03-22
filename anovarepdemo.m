%ANOVAREPDEMO
%This is a demo on the use of ANOVAREP function
%Example
%Run anovarepdemo
%
%           Created by Giuseppe Cardillo
%           giuseppe.cardillo-edta@poste.it
%
% To cite this file, this would be an appropriate format:
% Cardillo G. (2008) Anovarep: compute the Anova for repeated measures and
% Holm-Sidak test for multiple comparisons if Anova is positive. 
% http://www.mathworks.com/matlabcentral/fileexchange/18746

clear 
clc
load anrepdemo
disp('Biol Psychiatry, 2000; 47:89-1901')
disp('Testosterone levels (ng/dL)')
disp(array2table(x,'VariableNames',{'Soldier' 'At_begin_of_training' 'At_seize' 'h_12_after_seize' 'h_48_after_seize'}))
anovarep(x(2:end,:))
