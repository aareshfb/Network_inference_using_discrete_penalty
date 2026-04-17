%This is code was written by Aaresh to evaluate L1 algorithm with L0
%I have added the calculation for BIC

p=65;
n;
mask_diag= true;
T=3;
tol = 1e-6;



%% --- Paths (match Python defaults)
base_path = "../";                     % where *.table live
grn_path  = "../Simulated_GRNs/";      % where regulators idx + precision live

tol = 1e-6; % same tolerance style as your MATLAB BIC block

%% ---- Load regulator indices (same as Python regulator_idx) ----
idx = regulator_idx_matlab(grn_path);  % struct with fields hsc/cmp/gmp
reg_idx = idx.hsc;                     % regulators are same
%mu_list = [1e-6,1e-5,1e-4,1e-3,1e-2]; 
mu_list = [1e-4]
thresh = logspace(0,-8, 200); %Sparsity Coefficent
gamma = [1]; %Similarity Coefficent

data=load_scMTNI(n,base_path);

%% ---- Cluster similarity matrix ----
Adj = eye(T) + diag(ones(T-1,1), +1);    % 3x3 with superdiagonal ones

% ---- Load true precision matrices (same as Python read_prec) ----
true_networks = read_scMTNI_prec(grn_path, T);  % cell array {1..t}, each p x p

output_dir = 'Output/AUPRC';  % build folder path
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
n_str = strrep(num2str(n), '.', 'd');
fname = fullfile(output_dir, sprintf('l1_AUPRC_n%s', n_str));

if isfile(strcat(fname,'.csv'))
    fid = fopen([fname, '.csv'], 'at');
else 
    fid = fopen([fname, '.csv'], 'at');
    fprintf(fid, 'n,p,lambda(sparsity),mu0(thresh),gamma,bic,recall1,recall2,recall3,recall_avg,precision1,precision2,precision3,precision_avg,f1_1,f1_2,f1_3,f1_avg,Time(mins), mask_diag, no_edges, \n');
end

for gg = 1:length(gamma)
    gamma0=gamma(gg);
    for ii = 1:length(mu_list)


        mu0 = mu_list(ii)

        for l=1:T
            N(l) = size(data{l},2); % number of samples for each cluster
            mu(l) = mu0*sqrt(log(p)/N(l))*norm(data{l}); % the threshold for obtaining the backward mapping. Theoretically, it should scale with sqrt(log(p)/N).
            %The extra term is normalize the threshold level with the norm of sample covariance matrices
        end

        [B] = genBackward(data, mu); % B{l} is equal to the inverse of the soft thresholded sample covariance matrix, acting as the backward mapping    

        fprintf('number of edges after soft-thresholding \n')
        %for l = 1:L    
        %    fprintf(fid, '%d \n',  nnz(B{l}) - p);
        %end

        for i = 1:length(thresh)


            lambda0 = thresh(i)
            fprintf('parameters: mu0 : %d lambda0 : %d \n',mu0,lambda0);

            lambda = zeros(T,1);
            for l=1:T
                lambda(l) = lambda0*sqrt(N(l)*log(p)); % lambda(l) is the coefficient of ||theta^{(l)}||_1 in the objective. We set these values proportional to sqrt(N*log(p))
            end


            %set l1 and l2 penalty for network similarity
            Adj2 = gamma0*Adj;

            disp('Starting the estimation...')
            tic
            [Theta_l1] = elem_l1_parallel(B, Adj, lambda, N);
            Time = toc;       
            fprintf('Estimation done, time = %d mins \n', Time/60)


            %%number of non-zero edges in network
            no_edges=0
            for l = 1:T
                x = Theta_l1{l};
                no_edges=no_edges+nnz(abs(x)>1e-8) - p
            end
            no_edges=no_edges/2
            
            %% apply mask (remove diagonal)
            for i = 1:T
                mask = ~eye(p);
                true_networks2{i} = true_networks{i} .* mask;
                Theta_l1_2{i}     = Theta_l1{i} .* mask;
            end
            
            %% collect error metrics
            Errors = zeros(T,5);
            for i = 1:T
                [recall, precision, f1, max_error, norm_error] = error_cal(Theta_l1_2{i}(reg_idx,:), true_networks2{i}(reg_idx,:));
                Errors(i,:) = [recall, precision, f1, max_error, norm_error];
            end
            result = mean(Errors);
            disp(result)
            recalls = Errors(:,1);
            precisions = Errors(:,2);
            f1s = Errors(:,3);


            %% bic
            Theta_SPD = {};
            fprintf('finding nearest PSD matrix \n', Time)
            for t = 1:T
                try Theta_SPD{t} = nearestPSD(Theta_l1{t});
                    disp('Calculating nearestPSD')
                catch ME
                    Theta_SPD{t} = nearestPSD(full(Theta_l1{t}));
                end
            end
            
            bic = 0;
            total_nnz = 0;
            for t = 1:T                
                temp = Theta_SPD{t};
                temp(abs(temp) < tol) = 0;
                b = nnz(triu(temp,1));
                total_nnz = total_nnz + b;
                bic = bic + N(t)*(-log(det(Theta_SPD{t})) + ...
                    trace(cov(data{t}') * Theta_SPD{t}))+b*log(N(t))+4*b*log(p);            
            end
            
            %% Save the results to file.
            fprintf(fid,'%d,%d,%e,%e,%e,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d \n', ...
                n,p,lambda0,mu0,gamma0,bic, ...
                recalls(1),recalls(2),recalls(3),mean(recalls), ...
                precisions(1),precisions(2),precisions(3),mean(precisions), ...
                f1s(1),f1s(2),f1s(3),mean(f1s), ...
                Time/60,mask_diag,no_edges);

        end

    end
end
fclose(fid);





