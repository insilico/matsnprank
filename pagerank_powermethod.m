function r = pagerank_powermethod(names,G,p,resultsbase,capturedata)
% PAGERANK  Google's PageRank
% [header100,data100]=parsefile('Titers.Matrix.top100.txt');
% pagerank_powermethod(header100,data100,.85,'results');

% Uses URLs (SNPs) U and adjacency matrix G produced by parsefile,
% together with a damping factory p, (default is .85), to compute and plot
% a bar graph of page rank, and print the dominant URLs (SNPs) in page rank order.
% x = pagerank(U,G,p) returns the page ranks instead of printing.
% See also SURFER, SPY.

% Set defaults for optional params
if nargin < 5
    p = .85; 
    resultsbase = 'results';
    capturedata = false;
end

% G diagonal is information gain.  
Gdiag = diag(G);
% Sum of diagonal elements
Gtrace = sum(Gdiag);
% Creates a new matrix of 0s and sets Gs diagonal as the diagonal, 
% and subtracts from G.  Thus the result is G minus its diagonal.
Gnodiag = G - diag(Gdiag);   % sometimes not necessary

% colsum = out-degree, rowsum = in-degree
[n,n] = size(G);
colsum = sum(G, 1);  % 1xn vector of column sums
rowsum = sum(G, 2);  % nx1 vector of row sums

% Scale G column sums to be 1 (or 0 where there are no out links).
colsum_idx = find(colsum ~= 0);  % indices of colsum vector that are not zero 

% nxn sparse matrix with values 1/colsum for nonzero elements of colsum (indices
% given by k).  Other elements are zero.
D = sparse(colsum_idx, colsum_idx, Gtrace ./ colsum(colsum_idx), n, n);  

% compute initial T, Markov chain transition matrix
z_T(colsum_idx) = (1-p);
T = (p * G * D + Gdiag * z_T) / Gtrace;


% iterate power method
unit = ones(n,1); % column vector of 1s
r = unit/n; % initial arbitrary vector
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
    saveas(h, [resultsbase '-bar' num2str(p) '.eps'], 'psc2');

    % Scatter plot of SNPRank score vs. node degree
    scorevdeg = scatter(rowsum, r);
    xlabel('degree');
    ylabel('SNPRank score');
    title([resultsbase ' score vs. degree, p = ' num2str(p)])
    saveas(scorevdeg, [resultsbase '-degree-scatter' num2str(p) '.eps'], 'psc2')

    % Scatter plot of SNPRank score vs. information gain (IG)
    scorevig = scatter(Gdiag, r);
    xlabel('Information gain (IG)');
    ylabel('SNPRank score');
    title([resultsbase ' score vs. IG, p = ' num2str(p)])
    saveas(scorevig, [resultsbase '-IG-scatter' num2str(p) '.eps'], 'psc2')


    % Print SNPs in SNPrank order.
    fid = fopen([resultsbase '-' num2str(p) '.txt'], 'w');
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