function idx = regulator_idx_matlab(grn_path)
    % Matches Python regulator_idx() reading integer indices.
    % Python likely stores 0-based indices -> convert to 1-based for MATLAB.

    cell_types = { "hsc", "cmp", "gmp" };
    idx = struct();

    for k = 1:numel(cell_types)
        ct = cell_types{k};
        infile = grn_path + ct + "_regulators_idx.txt";
        v = readmatrix(infile, "FileType","text");
        v = v(:);

        % Convert 0-based -> 1-based (if file is already 1-based, this will be wrong)
        v = v + 1;

        idx.(ct) = v;
    end
end