function [theta_l1] = elem_l1_parallel(B, Adj, lambda, N)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Our goal is to estimate a set of sparse inverse covariance matrices for different clusters,
% whose similarities are controlled by an adjacency matrix: if two inverse
% covariance matrices are very similar, they differ only in a few elements.


% for each element (i,j) of the unknown inverse covariance matrix, we solve
% the following quadratic program:

% min. \sum_l (theta_{ij}^{(l)} - B_{ij}^{(l)})^2 + \sum_{k,l} W_{lk}|theta_{ij}^{(l)}-theta_{ij}^{(k)}| + \sum_l lambda |theta_{ij}^{(l)}| 

% this can be reformulated as a quadratic program with linear constraints, and can be solved efficiently using the "quadprog" function in MATLAB. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    p = size(B{1},1);
    L = length(B);

    [idx, idy, w] = find(Adj); % vectorizing the weights
    k = length(w);
    %% constructing the matrix for the quadratic term in the objective
    H = sparse(L+k+L,L+k+L);
    H(1:L,1:L) = sparse(2*diag(ones(size(N))));

    %% Constructing the linear inequality constraint Ax <= b
    A = sparse(2*k+2*L, L+k+L);
    b = sparse(size(A,1),1);

    A(sub2ind(size(A), 1:k, idx')) = 1;
    A(sub2ind(size(A), 1:k, idy')) = -1;
    A(sub2ind(size(A), 1:k, L+(1:k))) = -1;

    A(sub2ind(size(A), k + (1:k), idx')) = -1;
    A(sub2ind(size(A), k + (1:k), idy')) = 1;
    A(sub2ind(size(A), k + (1:k), L+(1:k))) = -1;

    A(sub2ind(size(A), 2*k+(1:L), 1:L)) = 1;
    A(sub2ind(size(A), 2*k+(1:L), L+k+(1:L))) = -1;

    A(sub2ind(size(A), 2*k+L+(1:L), 1:L)) = -1;
    A(sub2ind(size(A), 2*k+L+(1:L), L+k+(1:L))) = -1;

%%
    theta_l1 = cell(L,1);
    for l = 1:L
        theta_l1{l} = sparse(p,p);
    end
    temp_theta_l1 = [];
%     %%%%%%%%%%
    B = reshape(full(cell2mat(B)),p,p,L);
    F = zeros(L+length(w)+length(lambda),p*(p+1)/2);
    count = 1;
    for ii= 1:p
        for jj= ii:p
            f = squeeze(B(ii,jj,:));
            %%%%%%%%%%%%%%%
            if ii==jj
                lambda0 = zeros(size(lambda));
                F(:,count) = [-2*f; w; lambda0];
            else
                F(:,count) = [-2*f; w; lambda];
            end
%             F(:,count) = [f; w; lambda];
            %%%%%%%%%%%%%%%%%
            count = count + 1;
        end
    end
    %%%%%%%%%%
    parfor k = 1:size(F,2)
        if mod(k,500) == 0
            fprintf('solving for variable %d\n', k);
        end
        f = F(:,k);
        options = optimoptions('quadprog','Display','off');
        x = quadprog(H,f,A,b, [],[],[],[],[], options);
        temp_theta_l1 = [temp_theta_l1; x(1:L)'];
    end
    
    for l = 1:L
        temp = zeros(p);
        temp(tril(ones(p),0)==1) = temp_theta_l1(:,l);
        theta_l1{l} = sparse(temp)';
    end
    

    for l = 1:L
        theta_l1{l}(abs(theta_l1{l})<=1e-5) = 0; 
        theta_l1{l} = triu(theta_l1{l},1)+triu(theta_l1{l})';
    end

end



 