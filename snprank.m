function r = snprank(datafile, gamma, capturedata)
% SNPRank - SNP ranking algorithm

% Uses SNP names (SNPs) and adjacency matrix G produced by parsefile,
% together with a damping factor gamma, (default is .85), to compute 
% SNPRank scores.
% Returns the final SNPRank score vector r, in same order as original 
% data matrix/SNP names
%
% Usage:  pagerank_powermethod('gain-matrix.txt');
% Authors:  Brett McKinney and Nick Davis
% Email:  brett.mckinney@gmail.com, nick@nickdavis.name

% Set defaults for optional params
if nargin < 3
    gamma = .85;
    capturedata = false;
end

% Use file prefix (everything preceding .ext) for saved data files
namesplit = regexp(datafile, '\....$', 'split', 'stringanchors');
resultsbase = char(namesplit(1));

% Parse data matrix file, store headers (SNPs) and data as vars
[SNPs, G] = parsefile(datafile);

% G diagonal is information gain  
Gdiag = diag(G);
% Sum of diagonal elements
Gtrace = sum(Gdiag);

% colsum = out-degree, rowsum = in-degree (in undirected graphs
% out-degree = in-degree)
[n,n] = size(G);

% 1 x n (row) vector of column sums (d_j in Eqs. 4, 5 from SNPRank paper)
colsum = sum(G, 1);  
% n x 1 (column) vector of row sums
rowsum = sum(G, 2);  

% indices of colsum vector that are non-zero
colsum_nzidx = find(colsum ~= 0);  

% n x n sparse matrix with values 1/colsum for nonzero elements of colsum 
% indices given by colsum_nzidx).  Other elements are zero.
D = sparse(colsum_nzidx, colsum_nzidx, 1 ./ colsum(colsum_nzidx), n, n);  

% initialize second term multiplier as a column vector of 1s
T_nz = ones(1, n);

% non-zero elements of colsum/d_j have (1 - gamma) in the numerator of the 
% second term (Eq. 5 from SNPRank paper)
T_nz(colsum_nzidx) = 1 - gamma;

% compute initial T, Markov chain transition matrix
% first term (gamma * G * D) is zero when d_j = 0
T = (gamma * G * D) + (Gdiag * T_nz) / Gtrace;

% iterate power method
unit = ones(n, 1); % column vector of 1s
r = unit / n; % initial arbitrary vector
for i=1:5
    r = T * r;
    lambda = sum(r);  % temporary eigenvalue
    r = r / lambda;     % Normalize eivenvector so that sum(r) == 1.
end

% Bar graph of SNPRank, with labels for top SNPs
figure(1)
h = bar(r);
title('SNPRank scores')
% Capture top 10 SNPs by SNPRank, label vertically on bar graph
[sortedranks, snpindices] = sort(r, 'descend');
topsnps = sortedranks(1:10);
topindices = snpindices(1:10);
for i=1:length(topindices)
    % need to replace _ with \_ so Matlab doesn't subscript the labels
    text(topindices(i), topsnps(i),strrep(SNPs(topindices(i)),'_', '\_'), 'Rotation', 90)
end
% Only save figures and data to file if capturedata is set
if capturedata
    saveas(h, [resultsbase '-bar' num2str(gamma) '.eps'], 'psc2');
end

% Scatter plot of SNPRank score vs. node degree
figure(2)
scorevdeg = scatter(rowsum, r);
xlabel('degree');
ylabel('SNPRank score');
title([strrep(resultsbase,'_', '\_') ' score vs. degree, gamma = ' num2str(gamma)])
if capturedata
    saveas(scorevdeg, [resultsbase '-degree-scatter' num2str(gamma) '.eps'], 'psc2')
end

% Scatter plot of SNPRank score vs. information gain (IG)
figure(3)
scorevig = scatter(Gdiag, r);
xlabel('Information gain (IG)');
ylabel('SNPRank score');
title([strrep(resultsbase,'_', '\_') ' score vs. IG, gamma = ' num2str(gamma)])
if capturedata
    saveas(scorevig, [resultsbase '-IG-scatter' num2str(gamma) '.eps'], 'psc2')
end

% Print SNPs in SNPrank order
fid = 1; % unless capturedata set, only print to console
if capturedata
    fid = fopen([resultsbase '-' num2str(gamma) '.txt'], 'w');
end
[~,q] = sort(-r);
fprintf(fid,'SNP \t snp-rank \t IG (main effect) \t in \t out\n');
for k = 1:n
  j = q(k);
  fprintf(fid,'%s \t %8.4f \t %8.4f \t %4.0f \t %4.0f \t \n', ...
     SNPs{j}, r(j), G(j,j), rowsum(j), colsum(j));
end
if capturedata
    fclose(fid);
end
