function snprank(datafile, gamma, capturedata, showgraphs)
% SNPrank - SNP ranking algorithm

% Uses SNP names (SNPs) and adjacency matrix G, together with a 
% damping factor gamma, (default is .85), to compute SNPrank scores.
% Prints a series of rows containing the SNP name, SNPrank score,
% information gain, sorted in descending order by SNPrank.
%
% Usage:  snprank('gain-matrix.txt');
% Authors:  Brett McKinney, Nick Davis, and Ahwan Pandey
% Email:  brett.mckinney@gmail.com, nick@nickdavis.name, 
%         ahwan-pandey@utulsa.edu

% Set defaults for optional params, don't write results files or show
% graphs
if nargin < 4
    gamma = .85;
    capturedata = false;
    showgraphs = false;
end

% Use file prefix (everything preceding .ext) for saved data files
namesplit = regexp(datafile, '\....$', 'split', 'stringanchors');
resultsbase = char(namesplit(1));

% Parse data matrix file, store headers (SNPs) and data as vars
dataclass = importdata(datafile);
SNPs = dataclass.colheaders;
G = dataclass.data;

% G diagonal is information gain  
Gdiag = diag(G);
% Sum of diagonal elements
Gtrace = sum(Gdiag);

% dimensions of data matrix
[n,n] = size(G);

% colsum = degree of each SNP(in undirected graphs, out-degree = in-degree)
% 1 x n (row) vector of column sums (d_j in Eqs. 4, 5 from SNPRank paper)
colsum = sum(G, 1);  

% indices of colsum vector that are non-zero
colsum_nzidx = find(colsum ~= 0);  

% n x n sparse matrix with values 1/colsum for nonzero elements of colsum 
% indices given by colsum_nzidx).  Other elements are zero.
D = zeros(n);
for i = colsum_nzidx
	D(i,i) = 1 / colsum(colsum_nzidx(i));
end

% initialize second term multiplier as a row vector of 1s
T_nz = ones(1, n);

% non-zero elements of colsum/d_j have (1 - gamma) in the numerator of the 
% second term (Eq. 5 from SNPRank paper)
T_nz(colsum_nzidx) = 1 - gamma;

G_gpu = gsingle(G);
D_gpu = gsingle(D);

mult_in_gpu = G_gpu * D_gpu;

% compute initial T, Markov chain transition matrix
% first term (gamma * G * D) is zero when d_j = 0
T = (gamma * single(mult_in_gpu)) + (Gdiag * T_nz) / Gtrace;

% iterate power method
unit = ones(n, 1); % column vector of 1s
r = unit / n; % initial arbitrary vector
threshold = 10^-4; % cutoff for matrix convergence
converged = 0;
while ~converged
    r_old = r;
    r = T * r;
    lambda = sum(r);  % temporary eigenvalue
    r = r / lambda;     % Normalize eivenvector so that sum(r) == 1.
    if abs(r - r_old) < threshold
        converged = 1;
    end
end

if showgraphs
    % Bar graph of SNPrank, with labels for top SNPs
    figure(1)
    h = bar(r);
    title('SNPRank scores')
    % Capture top 10 SNPs by SNPrank, label vertically on bar graph
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

    % Scatter plot of SNPrank score vs. node degree
    figure(2)
    scorevdeg = scatter(colsum, r);
    xlabel('degree');
    ylabel('SNPRank score');
    title([strrep(resultsbase,'_', '\_') ' score vs. degree, gamma = ' num2str(gamma)])
    if capturedata
        saveas(scorevdeg, [resultsbase '-degree-scatter' num2str(gamma) '.eps'], 'psc2')
    end

    % Scatter plot of SNPrank score vs. information gain (IG)
    figure(3)
    scorevig = scatter(Gdiag, r);
    xlabel('Information gain (IG)');
    ylabel('SNPRank score');
    title([strrep(resultsbase,'_', '\_') ' score vs. IG, gamma = ' num2str(gamma)])
    if capturedata
        saveas(scorevig, [resultsbase '-IG-scatter' num2str(gamma) '.eps'], 'psc2')
    end
end

% Print SNPs in SNPrank order
fid = 1; % unless capturedata set, only print to console
if capturedata
    fid = fopen([resultsbase '-' num2str(gamma) '.txt'], 'w');
end
[~,q] = sort(-r);
fprintf(fid,'SNP\tSNPrank\tIG\n');
for k = 1:n
  j = q(k);
  fprintf(fid,'%s\t%8.6f\t%8.6f\n', ...
     SNPs{j}, r(j), G(j,j));
end
if capturedata
    fclose(fid);
end
