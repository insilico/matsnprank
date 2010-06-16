% Calculates a p-value for SNPRank scores
% Author:  Nick Davis
function calc_snprank_pval(ecfilename, matfilebase, numpermutes, resultbase)

% calculate max snprank score for EC data
[echdr, ecdata] = parsefile(ecfilename);
ec_snprank = pagerank_powermethod(echdr, ecdata, 0.85, resultbase, false);
ec_snprank
obs_max = max(ec_snprank);
max_perm_snprank = zeros(numpermutes,1);
for i=1:numpermutes
    [hdr, data] = parsefile([matfilebase '-' num2str(i) '.tab.txt']);
    snprank = pagerank_powermethod(hdr, data, 0.85, resultbase, false);
    max_perm_snprank(i) = max(snprank);
end
max_perm_snprank
pvalcalc = sum(max_perm_snprank >= obs_max) / numpermutes;
fprintf(1, 'Max EC SNPRank: %f\n', obs_max);
fprintf('Calculated p-val: %f\n', pvalcalc); 
