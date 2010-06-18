function r = pagerank_powermethod(names, G, gamma, resultsbase, capturedata)
% SNPRank - SNP ranking algorithm
% [header100,data100] = parsefile('Titers.Matrix.top100.txt');
% pagerank_powermethod(header100, data100);

% Uses names (SNPs) and adjacency matrix G produced by parsefile,
% together with a damping factory gamma, (default is .85), to compute and plot
% a bar graph of page rank, and print the dominant SNPs in rank order.
% Returns the final SNPRank score vector r

% Set defaults for optional params
if nargin < 5
    gamma = .85; 
    resultsbase = 'results';
    capturedata = false;
end

% G diagonal is information gain  
Gdiag = diag(G);
% Sum of diagonal elements
Gtrace = sum(Gdiag);

% colsum = out-degree, rowsum = in-degree (in undirected graphs
% out-degrees = in-degree)
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

% Only show figures and save data to file if capturedata is set
if capturedata
    % Bar graph of page rank, with labels for top SNPs
    shg
    h = bar(r);
    title('SNP Rank')

    [sortedranks, snpindices] = sort(r, 'descend');
    % Capture top 10 SNPs by SNPRank
    topsnps = sortedranks(1:10);
    topindices = snpindices(1:10);
    for i=1:length(topindices)
        % need to replace _ with \_ so Matlab doesn't subscript the labels
        text(topindices(i), topsnps(i),strrep(names(topindices(i)),'_', '\_'), 'Rotation', 90)
    end
    saveas(h, [resultsbase '-bar' num2str(gamma) '.eps'], 'psc2');

    % Scatter plot of SNPRank score vs. node degree
    scorevdeg = scatter(rowsum, r);
    xlabel('degree');
    ylabel('SNPRank score');
    title([resultsbase ' score vs. degree, gamma = ' num2str(gamma)])
    saveas(scorevdeg, [resultsbase '-degree-scatter' num2str(gamma) '.eps'], 'psc2')

    % Scatter plot of SNPRank score vs. information gain (IG)
    scorevig = scatter(Gdiag, r);
    xlabel('Information gain (IG)');
    ylabel('SNPRank score');
    title([resultsbase ' score vs. IG, gamma = ' num2str(gamma)])
    saveas(scorevig, [resultsbase '-IG-scatter' num2str(gamma) '.eps'], 'psc2')


    % Print SNPs in SNPrank order.
    fid = fopen([resultsbase '-' num2str(gamma) '.txt'], 'w');
    if nargout < 1
       [ignore,q] = sort(-r);
       fprintf(fid,'SNP \t snp-rank \t in \t out\n');
       fprintf('index \t snp-rank \t in \t out \t SNP \t IG (main effect)')
       k = 1;
       fprintf('n: %i', n)
       while (k <= n)
          j = q(k);
          fprintf(' %3.0f \t %8.4f \t %4.0f \t %4.0f \t %s \t %8.4f\n', ...
             j,r(j),rowsum(j),colsum(j),names{j},G(j,j))
         fprintf(fid,'%s \t %8.4f \t %4.0f \t %4.0f \t \n', ...
             names{j},r(j),rowsum(j),colsum(j));
          k = k + 1;
       end
    end
    fclose(fid);
end