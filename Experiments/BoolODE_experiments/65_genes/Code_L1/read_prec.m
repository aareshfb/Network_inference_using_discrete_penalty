function prec = read_prec()

base_path = fullfile('..', 'Ground_truth');
prec = cell(8,1);

expected_regs = arrayfun(@(i) sprintf('reg%d', i), 1:15, 'UniformOutput', false);
expected_genes = arrayfun(@(i) sprintf('gene%d', i), 1:50, 'UniformOutput', false);
ordered_cols = [expected_regs, expected_genes];

for i = 1:8
    filename = sprintf('network-%d.csv', i);
    path = fullfile(base_path, filename);

    T = readtable(path, 'ReadRowNames', true, 'VariableNamingRule', 'preserve');

    row_names = T.Properties.RowNames;
    col_names = T.Properties.VariableNames;

    if ~isequal(row_names, expected_regs')
        error('%s: row order mismatch.\nExpected rows: %s\nFound rows: %s', ...
            path, strjoin(expected_regs, ', '), strjoin(row_names', ', '));
    end

    if ~isequal(col_names, ordered_cols)
        error('%s: column order mismatch.\nExpected cols: %s\nFound cols: %s', ...
            path, strjoin(ordered_cols, ', '), strjoin(col_names, ', '));
    end

    arr = table2array(T);

    if ~isequal(size(arr), [15, 65])
        error('%s: expected shape (15, 65), found (%d, %d)', ...
            path, size(arr,1), size(arr,2));
    end

    prec{i} = arr;
    fprintf('Loaded %s: shape=(%d,%d)\n', filename, size(arr,1), size(arr,2));
end

end