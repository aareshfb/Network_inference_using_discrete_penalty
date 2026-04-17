function S = threshCov(x, tau, Sparse)
% Let M = x*x'. S is the soft-thresholded version of M with threshold tau.
% Input: x   -- n x m matrix.
%        tau -- threshold

% If we allocate more than this number times n nonzeros, terminate
% prematurely.
        
MEMORY_LIMIT = 1000; % x n

nn = size(x,2);
% Get partition size
n = size(x,1);
x2 = x;
skip = ceil(sqrt(n*MEMORY_LIMIT));
% if n <= 10*MEMORY_LIMIT 
%     % If problem size small, just do it explicitly
%     M = triu(x*x');
%     [ii_pos, jj_pos, kk_pos] = find(max(M-tau,0));
%     [ii_neg, jj_neg, kk_neg] = find(min(M+tau,0));
%     S = sparse([ii_pos; ii_neg],[jj_pos; jj_neg],[kk_pos; kk_neg],...
%                n,n);
%     S = S + triu(S,1)';
%    return
% end

% Partition x into cells
tmp = floor(n/skip);
rows = [repmat(skip,tmp,1); n - tmp*skip];
cumrows = [0;cumsum(rows)];
x = mat2cell(x, rows);

flag = 0;

switch nargin
    case 3
        flag = 1;
        Sparse = mat2cell(Sparse, rows, rows);
end

p = numel(x);

% Threshold one submatrix of M = x*x' at a time.
[ii,jj,kk] = deal(cell(1,p)); 
total_nz = 0;
for j = 1:p
    [this_ii, this_jj, this_kk] = deal(cell(1,j));
    for i = 1:j
        % Form matrix
        Mij = x{i}*x{j}'/nn;
        if i == j
            Mij = triu(Mij);
        end
        % Do soft-thresholding
        if flag == 1
            [ii_pos, jj_pos, kk_pos] = find(max(Mij-tau,0).*Sparse{i,j});
            [ii_neg, jj_neg, kk_neg] = find(min(Mij+tau,0).*Sparse{i,j});
        else
            [ii_pos, jj_pos, kk_pos] = find(max(Mij-tau,0));
            [ii_neg, jj_neg, kk_neg] = find(min(Mij+tau,0));
        end
        % Record nonzeros
        this_ii{i} = [ii_pos; ii_neg] + cumrows(i);
        this_jj{i} = [jj_pos; jj_neg] + cumrows(j);
        this_kk{i} = [kk_pos; kk_neg];
        % Sum nonzeros
        total_nz = total_nz + numel(this_ii{i});
        % Check for memory issues
        if total_nz > 10*MEMORY_LIMIT * n
            error('REACHED MEMORY LIMIT. EXITING....');
        end
    end
    % Assemble this column
    ii{j} = cat(1, this_ii{:});
    jj{j} = cat(1, this_jj{:});
    kk{j} = cat(1, this_kk{:});
end
% Assemble all columns
ii = cat(1, ii{:});
jj = cat(1, jj{:});
kk = cat(1, kk{:});

% form sparse matrix
S = sparse(ii,jj,kk,n,n);
S = S + triu(S,1)';
S = S-diag(diag(S))+diag(sparse(sum(x2.^2,2)))/nn;
end