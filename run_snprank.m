% Automates processing of a range of values for SNPrank
% Author:  Nick Davis
function run_snprank(matfile, resultbase)
pvals = [0 .25 .5 .85 1];
[hdr, data] = parsefile(matfile);
for i=1:length(pvals)
   pagerank_powermethod(hdr, data, pvals(i), resultbase, true) 
end
