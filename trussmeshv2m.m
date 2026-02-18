function [node, elem, face] = trussmeshv2m(vol, opt)
% [node,elem,face]=trussmeshv2m(vol,opt)
%
% Wrapper for trussmesh mesher (standalone executable version)
% Convert a binary (or multi-valued) volume to tetrahedral mesh
%
% -- this function is part of iso2mesh toolbox (http://iso2mesh.sf.net)

fprintf(1, 'creating surface and tetrahedral mesh from a multi-domain volume ...\n');

if ~(islogical(vol) || isa(vol, 'uint8')), error('trussmesher can only handle uint8 volumes'); end
if ~any(vol(:)), error('no labeled regions found in the input volume'); end
if islogical(vol), vol = uint8(vol) * 255; end
if max(vol(:)) > 63, error('trussmeshv2m supports at most 64 labels (0-63)'); end

% Parse options with defaults
defaults = struct('vertexspacing', 4, 'trussfactor', 1.3, 'springconstant', 0.1, ...
                  'compressionfactor', 1.3, 'terminatecriterion', 0.001, 'maxiter', 100000, ...
                  'debug', 1, 'surf_tet_qual_thresh', 0.6, 'retriangulateiterations', 50);
if nargin < 2, opt = struct(); end
if ~isstruct(opt), opt = struct('vertexspacing', opt); end
for f = fieldnames(defaults)', if ~isfield(opt, f{1}), opt.(f{1}) = defaults.(f{1}); end, end

% Setup paths and write input files
inputfile = mwpath('pre_trussmesh.raw');
headerfile = mwpath('pre_trussmesh.hdr');
outputfile = mwpath('post_trussmesh.raw');
deletemeshfile(outputfile);

fid = fopen(inputfile, 'wb'); fwrite(fid, permute(vol, [3,2,1]), 'uint8'); fclose(fid);
fid = fopen(headerfile, 'w'); fprintf(fid, '%d %d %d\n', size(vol, 1), size(vol, 2), size(vol, 3)); fclose(fid);

% Run mesher
cmd = sprintf('"%s" --input "%s" --header "%s" --output "%s" --spacing %f --truss %f --spring %f --compression %f --epsilon %f --maxiterations %d --retriangulateiterations %d', ...
    'trussmesh_cli', inputfile, headerfile, outputfile, ...
    opt.vertexspacing, opt.trussfactor, opt.springconstant, opt.compressionfactor, opt.terminatecriterion, opt.maxiter, opt.retriangulateiterations);
fprintf(1, "Trussmesh Command: %s\n", cmd);

if isfield(opt, 'log_file') && ~isempty(opt.log_file)
    cmd = sprintf('%s 2>&1 | tee "%s"', cmd, opt.log_file);
end
fprintf(1, "Trussmesh Command: %s\n", cmd);

[status, cmdout] = system(cmd);
if status ~= 0, error('trussmesh failed:\n%s', cmdout); end
fprintf(1, '%s', cmdout);
% if ~exist(outputfile, 'file'), error('output file not found'); end

% Load results
node = load_trussmesh_binary(outputfile);
fprintf(1, 'trussmesh: %d nodes\n', size(node, 1));

node_labels = extract_node_labels(node, vol);

% Warn about nodes that ended up entirely in the background
bg_only = node_labels == uint64(1);
if any(bg_only)
    bg_idx = find(bg_only);
    fprintf(2, 'WARNING: %d nodes have only background label (label 0):\n', length(bg_idx));
    for ii = 1:length(bg_idx)
        idx = bg_idx(ii);
        fprintf(2, '  node %d: [%.4f, %.4f, %.4f]\n', idx, node(idx,1), node(idx,2), node(idx,3));
    end
end

% Boundary nodes have multiple labels (more than one bit set) or touch background (bit 0)
is_boundary = (bitand(node_labels, node_labels - 1) ~= 0) | bitand(node_labels, uint64(1)) ~= 0;
print_node_summary(node_labels, is_boundary);

if size(node, 2) == 3  % 3D processing
    % Delaunay triangulation and remove long-edge tets

    fprintf(1, 'Performing delaunay triangulation...\n');
    elem = delaunay(node);
    fprintf(1, 'Done!\n');
    fprintf(1, 'Removing tets with long vertices...\n');
    tet_edges = nchoosek(1:4, 2);
    dists = zeros(size(elem,1), size(tet_edges,1));
    for ii = 1:size(tet_edges,1)
        dists(:,ii) = vecnorm(node(elem(:,tet_edges(ii,1)),:) - node(elem(:,tet_edges(ii,2)),:), 2, 2);
    end
    elem(any(dists > 2*opt.vertexspacing, 2), :) = [];
    fprintf(1, 'Done!\n');

    % Compute element labels via OR of vertex labels (boundary nodes zeroed)
    fprintf(1, 'Labeling tets...\n');
    tet_labels = node_labels(elem);
    tet_labels(is_boundary(elem)) = 0;
    elem(:,5) = double(bitor(bitor(tet_labels(:,1), tet_labels(:,2)), ...
                             bitor(tet_labels(:,3), tet_labels(:,4))));
    fprintf(1, 'Done!\n');

    % Remove low-quality unlabeled surface tets
    fprintf(1, 'Removing low quality tets...\n');
    qual = meshquality(node, elem);
    elem(qual < opt.surf_tet_qual_thresh & elem(:,5) == 0, :) = [];
    fprintf(1, 'Done!\n');

    % Label unlabeled tets if all neighbors share the same label
    fprintf(1, 'Labeling unlabeled tets...\n');
    conn = neighborelem(elem(:,1:3), size(node, 1));
    unlabeled = find(elem(:,5) == 0);
    for i = 1:length(unlabeled)
        tidx = unlabeled(i);
        nb = unique([conn{elem(tidx,1:3)}]);
        nb(nb == tidx) = [];
        if ~isempty(nb)
            nb_labels = elem(nb,5);
            nb_labels = nb_labels(nb_labels > 0);
            if ~isempty(nb_labels) && all(nb_labels == nb_labels(1))
                elem(tidx,5) = nb_labels(1);
            end
        end
    end
    fprintf(1, 'Done!\n');

    % Split spanning tets (tets with multiple labels)
    fprintf(1, 'Splitting spanning tets...\n');
    [node, elem, is_boundary] = split_spanning_tets(node, elem, node_labels, is_boundary);
    fprintf(1, 'Done!\n');

    % Label tets with all boundary vertices using centroid sampling
    all_bnd_tets = find(all(is_boundary(elem(:,1:4)), 2) & elem(:,5) == 0);
    for i = 1:length(all_bnd_tets)
        tidx = all_bnd_tets(i);
        centroid = mean(node(elem(tidx,1:4), :), 1);
        lbl = sample_vol(centroid, vol);
        elem(tidx,5) = double(label2bit(lbl));
    end

    % Verify tet orientations (ensure positive volume)
    fprintf(1, 'Verifying tet orientation...\n');
    for i = 1:size(elem,1)
        v = elem(i,1:4);
        tetvol = det([node(v(2),:)-node(v(1),:); node(v(3),:)-node(v(1),:); node(v(4),:)-node(v(1),:)]) / 6;
        if tetvol < 0, elem(i,1:4) = elem(i,[1 2 4 3]); end
    end
    fprintf(1, 'Done!\n');

    % Optional kite removal post-processing (disabled)
    % [node, elem] = remove_kites(node, elem);

    fprintf(1, 'Extracting faces...\n');
    % Extract faces per label, combining labels for shared faces
    labels = unique(elem(:,5));
    labels = labels(labels ~= 0);
    
    all_faces = [];
    all_labels = [];
    for ii = 1:length(labels)
        new_face = volface(elem(elem(:,5) == labels(ii), 1:4));
        all_faces = [all_faces; new_face];
        all_labels = [all_labels; uint64(labels(ii)) * ones(size(new_face,1), 1, 'uint64')];
    end
    
    % Sort vertices for grouping, then aggregate labels via bitor
    [face_sorted, sort_idx] = sortrows(sort(all_faces, 2));
    all_labels = all_labels(sort_idx);
    all_faces = all_faces(sort_idx, :);
    
    % Find unique faces and OR their labels together
    [~, ia, ic] = unique(face_sorted, 'rows', 'stable');
    face = all_faces(ia, :);
    face_labels = zeros(size(face,1), 1, 'uint64');
    for ii = 1:length(ic)
        face_labels(ic(ii)) = bitor(face_labels(ic(ii)), all_labels(ii));
    end
    face(:,4) = face_labels;
    fprintf(1, 'Done!\n');
else
    % TODO: 2D post-processing not yet implemented
    error('2D post-processing not yet supported');
end

fprintf(1, 'nodes: %d, triangles: %d, tetrahedra: %d\n', size(node, 1), size(face, 1), size(elem, 1));
end

%% Helper functions

function node = load_trussmesh_binary(filename)
    fid = fopen(filename, 'rb');
    nd = fread(fid, 1, 'uint32');
    npoints = fread(fid, 1, 'uint64');
    node = reshape(fread(fid, npoints * nd, 'float32'), [nd, npoints])';
    fclose(fid);
end

function print_node_summary(node_labels, is_boundary)
    % Print summary of internal and boundary nodes per label
    fprintf(1, '\n--- Node Summary ---\n');
    fprintf(1, 'Total nodes: %d\n', length(node_labels));
    fprintf(1, 'Total boundary: %d, Total internal: %d\n', sum(is_boundary), sum(~is_boundary));
    fprintf(1, '\nPer-label breakdown:\n');
    fprintf(1, '  Label | Internal | Boundary\n');
    fprintf(1, '  ------+----------+---------\n');
    
    % Check each possible label (0-63)
    for lbl = 0:63
        bit = bitshift(uint64(1), lbl);
        has_label = bitand(node_labels, bit) ~= 0;
        if ~any(has_label), continue; end
        
        n_internal = sum(has_label & ~is_boundary);
        n_boundary = sum(has_label & is_boundary);
        fprintf(1, '  %5d | %8d | %8d\n', lbl, n_internal, n_boundary);
    end
    fprintf(1, '--------------------\n\n');
end

function val = sample_vol(pt, vol)
    idx = min(max(round(pt) + 1, 1), size(vol));
    c = num2cell(idx(1:ndims(vol)));
    val = vol(c{:});
end

function mask = label2bit(lbl)
    % Convert label (0-63) to uint64 bitmask
    % Label 0 (background/outside) uses bit 0 (LSB)
    % Label N uses bit N
    if lbl >= 0 && lbl < 64
        mask = bitshift(uint64(1), lbl);
    else
        mask = uint64(0);
    end
end

function node_labels = extract_node_labels(node, vol)
    % Returns uint64 array where each bit represents a label
    % Bit 0 = background (label 0), Bit N = label N
    % Nodes with multiple bits set or bit 0 set are on boundaries
    ndim = size(node, 2);
    [offsets{1:ndim}] = ndgrid([-0.5, 0, 0.5]);
    offsets = reshape(cat(ndim+1, offsets{:}), [], ndim);
    
    % Vectorized: sample all nodes at all offsets, convert to bitmasks, OR together
    node_labels = zeros(size(node,1), 1, 'uint64');
    for j = 1:size(offsets, 1)
        idx = min(max(round(node + offsets(j,:)) + 1, 1), size(vol));
        sampled = vol(sub2ind(size(vol), idx(:,1), idx(:,2), idx(:,3)));
        node_labels = bitor(node_labels, bitshift(uint64(1), uint64(sampled)));
    end
end

function [node, elem, is_boundary] = split_spanning_tets(node, elem, node_labels, is_boundary)
    % Split tets that span multiple labels by inserting midpoint nodes
    spanning_tets = find(bitand(uint64(elem(:,5)), uint64(elem(:,5))-1) ~= 0);
    tet_labels = node_labels(elem(spanning_tets, 1:4));
    
    new_tets = [];
    new_nodes = [];
    tets_to_remove = [];
    
    for i = 1:length(spanning_tets)
        tidx = spanning_tets(i);
        tet_node_ids = elem(tidx, 1:4);
        tet_node_labels = tet_labels(i, :);
        
        is_bnd = is_boundary(tet_node_ids);
        bnd_nodes = tet_node_ids(is_bnd);
        int_nodes = tet_node_ids(~is_bnd);
        int_labels = tet_node_labels(~is_bnd);
        
        % Case: 2 boundary, 2 internal with different labels
        if sum(is_bnd) == 2 && length(int_nodes) == 2 && int_labels(1) ~= int_labels(2)
            mid_pt = (node(int_nodes(1),:) + node(int_nodes(2),:)) / 2;
            new_node_id = size(node,1) + size(new_nodes,1) + 1;
            new_nodes = [new_nodes; mid_pt];
            
            new_tets = [new_tets;
                bnd_nodes(1), bnd_nodes(2), int_nodes(1), new_node_id, int_labels(1);
                bnd_nodes(1), bnd_nodes(2), int_nodes(2), new_node_id, int_labels(2)];
            tets_to_remove = [tets_to_remove; tidx];
            
        % Case: 1 boundary, 3 internal spanning 2 labels
        elseif sum(is_bnd) == 1 && length(int_nodes) == 3
            unique_labels = unique(int_labels(int_labels > 0));
            if length(unique_labels) == 2
                [new_nodes, new_tets, tets_to_remove] = split_1bnd_3int(...
                    node, new_nodes, new_tets, tets_to_remove, tidx, ...
                    bnd_nodes(1), int_nodes, int_labels, unique_labels);
            end
            
        % Case: 0 boundary, 4 internal spanning 2 labels (2-2 split)
        elseif sum(is_bnd) == 0
            unique_labels = unique(int_labels(int_labels > 0));
            if length(unique_labels) == 2
                mask1 = int_labels == unique_labels(1);
                mask2 = int_labels == unique_labels(2);
                if sum(mask1) == 2 && sum(mask2) == 2
                    [new_nodes, new_tets, tets_to_remove] = split_0bnd_4int(...
                        node, new_nodes, new_tets, tets_to_remove, tidx, ...
                        int_nodes(mask1), int_nodes(mask2), unique_labels);
                end
            end
        end
    end
    
    % Update mesh
    if ~isempty(new_nodes)
        node = [node; new_nodes];
        is_boundary = [is_boundary; true(size(new_nodes,1), 1)];
    end
    elem(tets_to_remove, :) = [];
    if ~isempty(new_tets), elem = [elem; new_tets]; end
end

function [new_nodes, new_tets, tets_to_remove] = split_1bnd_3int(...
        node, new_nodes, new_tets, tets_to_remove, tidx, bnd_node, int_nodes, int_labels, unique_labels)
    % Split tet with 1 boundary node and 3 internal nodes spanning 2 labels
    mask1 = int_labels == unique_labels(1);
    mask2 = int_labels == unique_labels(2);
    nodes_l1 = int_nodes(mask1);
    nodes_l2 = int_nodes(mask2);
    
    % Create midpoints on cross-label edges
    edge_node_ids = [];
    for j = 1:length(nodes_l1)
        for k = 1:length(nodes_l2)
            mid_pt = (node(nodes_l1(j),:) + node(nodes_l2(k),:)) / 2;
            new_node_id = size(node,1) + size(new_nodes,1) + 1;
            new_nodes = [new_nodes; mid_pt];
            edge_node_ids = [edge_node_ids; new_node_id];
        end
    end
    
    if sum(mask1) == 2 && sum(mask2) == 1
        new_tets = [new_tets;
            bnd_node, nodes_l1(1), nodes_l1(2), edge_node_ids(1), unique_labels(1);
            bnd_node, nodes_l1(2), edge_node_ids(1), edge_node_ids(2), unique_labels(1);
            bnd_node, nodes_l2(1), edge_node_ids(1), edge_node_ids(2), unique_labels(2)];
    elseif sum(mask1) == 1 && sum(mask2) == 2
        new_tets = [new_tets;
            bnd_node, nodes_l1(1), edge_node_ids(1), edge_node_ids(2), unique_labels(1);
            bnd_node, nodes_l2(1), nodes_l2(2), edge_node_ids(1), unique_labels(2);
            bnd_node, nodes_l2(2), edge_node_ids(1), edge_node_ids(2), unique_labels(2)];
    end
    tets_to_remove = [tets_to_remove; tidx];
end

function [new_nodes, new_tets, tets_to_remove] = split_0bnd_4int(...
        node, new_nodes, new_tets, tets_to_remove, tidx, nodes_l1, nodes_l2, unique_labels)
    % Split internal tet with 2-2 label split using midpoints and centroid
    n1a = nodes_l1(1); n1b = nodes_l1(2);
    n2a = nodes_l2(1); n2b = nodes_l2(2);
    
    % Create midpoints on all 4 cross-label edges
    cross_edges = [n1a n2a; n1a n2b; n1b n2a; n1b n2b];
    mids = zeros(4,1);
    for j = 1:4
        mid_pt = (node(cross_edges(j,1),:) + node(cross_edges(j,2),:)) / 2;
        new_node_id = size(node,1) + size(new_nodes,1) + 1;
        new_nodes = [new_nodes; mid_pt];
        mids(j) = new_node_id;
    end
    
    % Create interface quad centroid
    quad_centroid = mean(new_nodes(end-3:end,:), 1);
    centroid_id = size(node,1) + size(new_nodes,1) + 1;
    new_nodes = [new_nodes; quad_centroid];
    
    % Split into tets for each label
    for j = 1:4
        new_tets = [new_tets; n1a, n1b, mids(j), centroid_id, unique_labels(1)];
        new_tets = [new_tets; n2a, n2b, mids(j), centroid_id, unique_labels(2)];
    end
    tets_to_remove = [tets_to_remove; tidx];
end

% function [node, elem] = remove_kites(node, elem)
% % Remove low-quality "kite" tetrahedra by inserting midpoint on longest diagonals
% qual_thresh = 0.5;
% fprintf("Tets: %d\n", size(elem,1));
% 
% qualities = meshquality(node, elem(:,1:4));
% % Ignore multi-label and unlabeled elements
% qualities(elem(:,5) == 0 | bitand(uint64(elem(:,5)), uint64(elem(:,5))-1)) = nan;
% neigh = neighbors(triangulation(elem(:,1:4), node));
% qualities(any(isnan(neigh), 2)) = nan;
% 
% while true
%     [m, i] = min(qualities);
%     if m > qual_thresh, break; end
%     qualities(i) = nan;  % Mark as processed
% 
%     tet = elem(i, 1:4);
%     fprintf("Quality: %f, index: %d\n", m, i);
%     
%     % Find edge lengths
%     tet_edges = nchoosek(1:4, 2);
%     dists = zeros(size(tet_edges,1), 1);
%     for ii = 1:size(tet_edges,1)
%         dists(ii) = norm(node(tet(tet_edges(ii,1)),:) - node(tet(tet_edges(ii,2)),:));
%     end
% 
%     % Get two longest edges (diagonals)
%     [~, ix] = sort(dists, 'descend');
%     diag0_loc = tet_edges(ix(1), :);
%     diag1_loc = tet_edges(ix(2), :);
%     diag0_glob = tet(diag0_loc);
%     diag1_glob = tet(diag1_loc);
% 
%     % Find shortest distance between two line segments using parametric representation
%     % Line 1: L1(s) = P1 + s*u, Line 2: L2(t) = P3 + t*v, where s,t ∈ [0,1]
%     P = [node(tet(diag0_loc), :); node(tet(diag1_loc), :)];
%     u = P(2,:) - P(1,:);  % Direction of first diagonal
%     v = P(4,:) - P(3,:);  % Direction of second diagonal
%     w = P(1,:) - P(3,:);
% 
%     a = dot(u,u); b = dot(u,v); c = dot(v,v); d = dot(u,w); e = dot(v,w);
%     denom = a*c - b*b;
%     
%     if abs(denom) < 1e-10  % Lines nearly parallel
%         s = 0; t = e/c;
%     else  % Lines are skew
%         s = (b*e - c*d) / denom;
%         t = (a*e - b*d) / denom;
%     end
%     s = max(0, min(1, s));
%     t = max(0, min(1, t));
%     
%     % Midpoint of shortest segment between diagonals
%     new_point = ((P(1,:) + s*u) + (P(3,:) + t*v)) / 2;
% 
%     % Find all tets containing either diagonal
%     e = elem(:, 1:4);
%     has_diag0 = any(e == diag0_glob(1), 2) & any(e == diag0_glob(2), 2);
%     has_diag1 = any(e == diag1_glob(1), 2) & any(e == diag1_glob(2), 2);
%     tets_to_split = find(has_diag0 | has_diag1);
% 
%     % Skip if tets span multiple labels (would mess with boundaries)
%     if numel(unique(elem(tets_to_split, 5))) ~= 1
%         fprintf("Skipping, not all tets are the same label\n");
%         continue;
%     end
% 
%     % Generate new tets by splitting
%     new_tets = [];
%     new_pt_idx = size(node, 1) + 1;
% 
%     for jj = 1:length(tets_to_split)
%         tet_idx = tets_to_split(jj);
%         curr_tet = elem(tet_idx, 1:4);
%         has_d0 = all(ismember(diag0_glob, curr_tet));
%         has_d1 = all(ismember(diag1_glob, curr_tet));
% 
%         if has_d0 && has_d1
%             continue;  % Original kite tet - remove without replacement
%         elseif has_d0
%             other = setdiff(curr_tet, diag0_glob);
%             new_tets = [new_tets; other(1), other(2), diag0_glob(1), new_pt_idx;
%                                   other(1), other(2), diag0_glob(2), new_pt_idx];
%         elseif has_d1
%             other = setdiff(curr_tet, diag1_glob);
%             new_tets = [new_tets; other(1), other(2), diag1_glob(1), new_pt_idx;
%                                   other(1), other(2), diag1_glob(2), new_pt_idx];
%         end
%     end
% 
%     % Check quality before committing
%     new_qual = meshquality([node; new_point], new_tets);
%     if any(new_qual < m) || any(isnan(new_qual))
%         fprintf("New tets have worse quality than original, skipping\n");
%         continue;
%     end
% 
%     % Commit changes
%     node = [node; new_point];
%     new_tets = [new_tets, elem(i,5) * ones(size(new_tets,1), 1)];
%     elem(tets_to_split, :) = [];
%     tets_to_split(tets_to_split > size(qualities,1)) = [];
%     qualities(tets_to_split) = [];
%     elem = [elem; new_tets];
% end
% end
