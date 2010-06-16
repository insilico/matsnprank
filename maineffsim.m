function maineffsim(base, ext, nummat, mu)
% Simulates main effect in an information gain matrix
% Base is the prefix of the file, ext is the extension,
% nummat is the number of randomly generated matrices
% mu is the value to use for the random deviates chosen from the normal dist
% E.g. for main1.mat, main2.mat, ..., main100.mat
% base = "main"
% ext = "mat"
% nummat = 100

% Determine size of matrices (assume all are the same dimensions)
[header, data] = parsefile(base + '1.mat');
[rows, cols] = size(data);
% preallocate A
A = zeros(rows, cols, nummat);

% read in A(i) matrices from file params passed
% run SNPRank on main effect matrices
for i = 1:nummat
    [header, data] = parsefile([base num2str(i) '.' + ext]);
    A(:,:,i) = data;
    % set 1,1 element of matrix to elevated main effect value
    A(1, 1, i) = mean(diag(A(i))) + normrnd(mu,1);
    pagerank_powermethod(header, A, mu);
end

