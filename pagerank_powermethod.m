function x = pagerank_powermethod(U,G,p,resultsbase,capturedata)
% PAGERANK  Google's PageRank
% [header100,data100]=parsefile('Titers.Matrix.top100.txt');
% pagerank_powermethod(header100,data100,.85,'results');

% Uses URLs (SNPs) U and adjacency matrix G produced by parsefile,
% together with a damping factory p, (default is .85), to compute and plot
% a bar graph of page rank, and print the dominant URLs (SNPs) in page rank order.
% x = pagerank(U,G,p) returns the page ranks instead of printing.
% See also SURFER, SPY.

if nargin < 3, p = .85; end

if nargin < 5, capturedata = true; end

% comment out because we want to keep the diagonal
% Eliminate any self-referential links
%Gdiag=diag(G)/sum(diag(G));
Gdiag=diag(G);
Gtrace=sum(Gdiag);
%Gdiag=Gdiag/Gtrace'
%G = G - diag(Gdiag);
Gtmp = G - diag(Gdiag);   % sometimes not necessary

% G diagonal is information gain.  
%Gij = InteractionInformation(SNP_i,SNP_j), i!=j

% examples
%G2=[2 1 0;0 1 2;0 2 1];
%U={'a','b','c'};

% c = out-degree, r = in-degree

[n,n] = size(G);
c = sum(G,1);  % 1xn vector of column sums
r = sum(G,2);  % nx1 vector of row sums

c2=sum(Gtmp,1);
% Scale column sums to be 1 (or 0 where there are no out links).

k = find(c~=0);    % indices of c vector that are not zero 
D = sparse(k,k,Gtrace./c(k),n,n);  % nxn matrix with values 1/c for
                              % nonzero elements of c (indices given by k)
                              % Other elements are zero.
% Solve (I - p*G*D)*x = e

%e = ones(n,1);    % original e
unit = ones(n,1);
e=unit.*Gdiag;      % modified to include IG
I = speye(n,n);
%z_T=ones(1,n)/n;    % original z_T
%z_T=e'/n;           % modfied to include IG
%z_T(k)=(1-p)/n;
z_T(k)=(1-p);
A=(p*G*D+e*z_T)/Gtrace;

% Instead of power method, this works
% x = (I - p*G*D)\e;
% x = Z\e is the solution to the equation Z*x = e computed by Gaussian
% elimination.

%find indices of A ~= 1
% A
% nononeindices = find(sum(A,1)~=1)
% for i=1:length(nononeindices)
%     A(:, nononeindices(i)) = A(:, nononeindices(i)) / sum(A(:, nononeindices(i)));
% end
% 
% A
% find(sum(A,1)~=1)
% A(sum(A,1)~=1)

% iterate power method
x=unit/n;
for i=1:5
    x=A*x;
    lambda=sum(x);  % temporary eigenvalue
    x=x/lambda;     % Normalize eivenvector so that sum(x) == 1.
end

lambda  % should be 1

% find blocks of graph
%[p,q,r,s,cc,rr]=dmperm(G);
%blocks=G(p,q)


% Only show figures and save data to file if capturedata is set
if capturedata
    % Bar graph of page rank, with labels for top SNPs
    shg
    h = bar(x);
    title('SNP Rank')

    [sortedranks, snpindices] = sort(x, 'descend');
    % Capture top 10 SNPs by SNPRank
    topsnps = sortedranks(1:10);
    topindices = snpindices(1:10);
    for i=1:length(topindices)
        % need to replace _ with \_ so Matlab doesn't subscript the labels
        text(topindices(i), topsnps(i),strrep(U(topindices(i)),'_', '\_'), 'Rotation', 90)
    end
    saveas(h, [resultsbase '-bar' num2str(p) '.eps'], 'psc2');

    % Scatter plot of SNPRank score vs. node degree
    scorevdeg = scatter(r,x);
    xlabel('degree');
    ylabel('SNPRank score');
    title([resultsbase ' score vs. degree, p = ' num2str(p)])
    saveas(scorevdeg, [resultsbase '-degree-scatter' num2str(p) '.eps'], 'psc2')

    % Scatter plot of SNPRank score vs. information gain (IG)
    scorevig = scatter(Gdiag,x);
    xlabel('Information gain (IG)');
    ylabel('SNPRank score');
    title([resultsbase ' score vs. IG, p = ' num2str(p)])
    saveas(scorevig, [resultsbase '-IG-scatter' num2str(p) '.eps'], 'psc2')


    % Print SNPs in SNPrank order.
    fid = fopen([resultsbase '-' num2str(p) '.txt'], 'w');
    if nargout < 1
       [ignore,q] = sort(-x);
       fprintf(fid,'S0NP \t snp-rank \t in \t out\n')
       disp(sprintf('index \t snp-rank \t in \t out \t SNP \t IG (main effect)'))
       k = 1;
       while (k <= n) %&& (x(q(k)) <= .005)
          j = q(k);
          disp(sprintf(' %3.0f \t %8.4f \t %4.0f \t %4.0f \t %s \t %8.4f', ...
             j,x(j),r(j),c(j),U{j},G(j,j)))
         %fprintf(fid,'%3.0f \t %8.4f \t %4.0f \t %4.0f \t %s \n', ...
         %    j,x(j),r(j),c(j),U{j});
         fprintf(fid,'%s \t %8.4f \t %4.0f \t %4.0f \t \n', ...
             U{j},x(j),r(j),c(j));
          k = k+1;
       end
    end
    fclose(fid);
end