function [recall, precision, f1, max_error, norm_error] = error_cal(A, invSigma)

    % ---- Sparsity pattern ----
    Sparse_estimation = sparse(abs(A) > 1e-8);
    Sparse_true = sparse(abs(invSigma) > 1e-8);

    % ---- Precision / Recall ----
    TP = sum(sum(Sparse_true .* Sparse_estimation));
    FP = sum(sum(Sparse_estimation)) - TP;
    FN = sum(sum(Sparse_true)) - TP;

    recall = TP / (TP + FN);
    precision = TP / (TP + FP);
    f1 = 2 * recall * precision / (recall + precision);

    % ---- Error metrics ----
    E = A - invSigma;

    [~,~,nnz_E] = find(E);
    [~,~,nnz_invSigma] = find(invSigma);

    max_error = max(abs(nnz_E));
    norm_error = norm(nnz_E) / norm(nnz_invSigma);

end