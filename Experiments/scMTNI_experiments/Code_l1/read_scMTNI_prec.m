function true_networks = read_scMTNI_prec(grn_path, t)
    % Matches Python read_prec: drop first row and first col from CSV.
    names = { "hsc", "cmp", "gmp" };
    true_networks = cell(1,t);

    for i = 1:t
        path = grn_path + names{i} + "_precision_matrix.csv";
        M = readmatrix(path,"NumHeaderLines",1);
        true_networks{i} = M(1:end, 2:end);
    end
end