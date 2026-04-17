function data = load_data()

    base_path = fullfile('..', 'Expression_data');

    % expected names
    expected_regs = arrayfun(@(i) sprintf('reg%d', i), 1:15, 'UniformOutput', false);
    expected_genes = arrayfun(@(i) sprintf('gene%d', i), 1:50, 'UniformOutput', false);
    ordered_rows = [expected_regs, expected_genes];

    data = cell(1,8);

    rng(0); % reproducibility

    for cluster_idx = 1:8

        filename = sprintf('cluster_%d_sparse_expression.csv', cluster_idx);
        path = fullfile(base_path, filename);

        % Read table with row names preserved
        T = readtable(path, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');

        % Convert to numeric array
        X = table2array(T);
        row_names = T.Properties.RowNames;

        % Find missing rows
        [~, loc] = ismember(ordered_rows, row_names);
        missing_idx = find(loc == 0);

        % Add missing rows with small Gaussian noise
        for k = 1:length(missing_idx)
            row_name = ordered_rows{missing_idx(k)};
            noise = 1e-5 * randn(1, size(X,2));

            X = [X; noise];
            row_names{end+1,1} = row_name;
        end

        % Recompute ordering after adding rows
        [tf, order] = ismember(ordered_rows, row_names);
        if ~all(tf)
            error('Some required rows are still missing in %s.', filename);
        end
        X = X(order, :);

        % Store as features x cells = 65 x n
        data{cluster_idx} = X;

        fprintf('Loaded %s: added missing=%d, final shape=(%d,%d)\n', ...
            filename, length(missing_idx), size(X,1), size(X,2));
    end

end