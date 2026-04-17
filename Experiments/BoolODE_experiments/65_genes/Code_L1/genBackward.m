function [B] = genBackward(x,mu)

tic

% x = Samples_small

L = size(x,2);
cond = zeros(1, L);
for l = 1:L
    S{l} = threshCov(x{l},mu(l));  
    B{l} = inv(S{l});
    cond(l) = condest(S{l});
    fprintf('condition number for covariance matrix %d is : %d\n', l, condest(S{l}));    
end

time = toc;
fprintf('Generating approximate backward mappings done, time = %d\n', time)

end



