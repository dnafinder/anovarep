# anovarep
# Analysis of variance for repeated measures. 
This function executes the analysis of variance when subjects underwent
several treatments. This function is similar to ANOVA2 Matlab Function,
but there are three differences:
1) the output of ANOVA table;
2) the graphical plot of anova table;
3) if p-value<alpha this function executes the Holm-Sidak test for multiple
comparison test to highlight differences between treatments.

Syntax: 	ANOVAREP(X,ALPHA)
     
    Inputs:
          X - data matrix. 
          ALPHA - significance level (default = 0.05).
    Outputs:
          - Anova table.
          - Graphical plot of anova table.
          - Holm-Sidak test (eventually)

     Example: anovarepdemo

          Created by Giuseppe Cardillo
          giuseppe.cardillo-edta@poste.it

To cite this file, this would be an appropriate format:
Cardillo G. (2008) Anovarep: compute the Anova for repeated measures and
Holm-Sidak test for multiple comparisons if Anova is positive. 
http://www.mathworks.com/matlabcentral/fileexchange/18746
