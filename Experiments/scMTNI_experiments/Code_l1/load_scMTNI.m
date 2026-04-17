function data = load_scMTNI(n, base_path)
    % Matches Python load_dataset:
    % df = read_csv(tab, header=None)
    % arr = df.iloc[1:,1:].to_numpy()
    % data[idx] = arr.T
    %
    % We then transpose once more so MATLAB uses p x n.

    t = 3;

    if n == 2000
        files = { base_path+"hsc.table", base_path+"cmp.table", base_path+"gmp.table" };
    elseif n == 1000
        files = { base_path+"hsc_n1000.table", base_path+"cmp_n1000.table", base_path+"gmp_n1000.table" };
    elseif n == 200
        files = { base_path+"hsc_n200.table", base_path+"cmp_n200.table", base_path+"gmp_n200.table" };
    else
        error("n must be one of {200, 1000, 2000}");
    end

    data = cell(1,t);
    for i = 1:t
        M = readmatrix(files{i}, "FileType","text", "Delimiter","\t");
        A = M(1:end, 2:end);     % drop first col
        data{i} = A;             % A is genes x cells
        data{i} = data{i};       % genes x cells
        fprintf('Loaded %s: %dx%d (genes x cells)\n', files{i}, size(A,1), size(A,2));
    end

    p = size(data{1}, 1);

    % MATLAB pipeline expects p x N, so this is already correct:
    % genes x cells = p x n.
end