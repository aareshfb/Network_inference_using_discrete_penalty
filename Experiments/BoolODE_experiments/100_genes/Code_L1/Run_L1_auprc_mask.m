%This is code was written by Aaresh to evaluate L1 algorithm with L0
%I have added the calculation for BIC

p=100;
n=2500;
mask_diag= true;
T=10;
tol = 1e-6;
no_reg=25;

tol = 1e-6; % same tolerance style as your MATLAB BIC block

mu_list = [1e-4];%Thresholding coefficent
thresh = logspace(0,-8, 200); %Sparsity
gamma = [1e-4]; %Similarity Coefficent

data=load_data(no_reg,p,T);

%% ---- Cluster similarity matrix ----
Adj = eye(10);
edges = [1 5; 5 10; 4 10; 4 7; 4 8; 9 10; 3 9; 6 9; 2 6];
Adj(sub2ind(size(Adj), edges(:,1), edges(:,2))) = 1;

% ---- Load true precision matrices (same as Python read_prec) ----
true_networks = read_prec(no_reg,p,T);  % cell array {1..t}, each p x p

output_dir = 'Output/AUPRC_1e-4';  % build folder path
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end
n_str = strrep(num2str(n), '.', 'd');
fname = fullfile(output_dir, sprintf('Code_L1_AUPRC_n%s', n_str));

if isfile(strcat(fname,'.csv'))
    fid = fopen([fname, '.csv'], 'at');
else 
    fid = fopen([fname, '.csv'], 'at');
    fprintf(fid, 'n,p,lambda(thest),mu0(sparsity),gamma,bic,recall_global,precision_global,f1_global,Time(mins), mask_diag, no_edges, \n');
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
                no_edges=no_edges+nnz(abs(x)>1e-5) - p
            end
            no_edges=no_edges/2
            
            %% apply mask (remove diagonal)
            for i = 1:T
                mask = ~eye(p);
                true_networks2{i} = true_networks{i};
                idx = sub2ind(size(true_networks2{i}), 1:no_reg, 1:no_reg);
                true_networks2{i}(idx) = 0;
                Theta_l1_2{i}     = Theta_l1{i} .* mask;
            end
            
            %% collect error metrics
            Errors = zeros(T,5);
            TP = 0;
            FP = 0;
            FN = 0;
            for i = 1:T
                pred  = Theta_l1_2{i}(1:no_reg,:);
                truth = true_networks2{i}(:,:);
            
                pred_bin  = abs(pred) > 1e-5;
                truth_bin = abs(truth) > 0;
            
                TP = TP + sum(sum(pred_bin & truth_bin));
                FP = FP + sum(sum(pred_bin & ~truth_bin));
                FN = FN + sum(sum(~pred_bin & truth_bin));
            
                % keep per-cluster metrics (unchanged)
                [recall, precision, f1, max_error, norm_error] = error_cal(pred, truth);
                Errors(i,:) = [recall, precision, f1, max_error, norm_error];
            end
            result = mean(Errors);
            disp(result)
            recalls = Errors(:,1);
            precisions = Errors(:,2);
            f1s = Errors(:,3);
            Precision_global = TP / (TP + FP + 1e-10);
            Recall_global    = TP / (TP + FN + 1e-10);
            F1_global        = 2 * Precision_global * Recall_global / (Precision_global + Recall_global + 1e-10);
            
            fprintf('GLOBAL Precision: %.4f\n', Precision_global);
            fprintf('GLOBAL Recall: %.4f\n', Recall_global);
            fprintf('GLOBAL F1: %.4f\n', F1_global);

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
            fprintf(fid,'%d,%d,%e,%e,%e,%.3f,%.3f,%.3f,%.3f,%.3f,%d,%d \n', ...
                n,p,lambda0,mu0,gamma0,bic, ...
                Recall_global, ...
                Precision_global, ...
                F1_global, ...
                Time/60,mask_diag,no_edges);

        end

    end
end
fclose(fid);





