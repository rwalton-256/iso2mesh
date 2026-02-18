function save_labeled_ply(filename, node, face, elem)
% SAVE_LABELED_PLY  Save a labeled surface mesh as a colored binary PLY file.
%
%   save_labeled_ply(filename, node, face, elem)
%
%   face(:,4) should contain label bitmasks (from trussmeshv2m).
%   elem(:,5) should contain element labels (used for tet PLY export).
%
%   This writes TWO files:
%     filename            - surface boundary mesh with per-face colors
%     [filename]_tet.ply  - tet mesh with per-cell region label as scalar
%
%   MeshLab can open the surface PLY directly. For the tet mesh,
%   consider Paraview (.vtu) instead — see save_labeled_vtu below.

%% --- Decode labels from bitmasks ---
% face(:,4) is a uint64 bitmask; convert to the dominant non-background label
% Use bitcmp for two's complement since MATLAB uint64 negation wraps wrong
face_bits = uint64(face(:,4));
shifted = bitshift(face_bits, -1);  % drop bit 0 (background)
face_labels = zeros(size(face,1), 1);
valid = shifted > 0;
neg_shifted = bitcmp(shifted(valid)) + uint64(1);
lowest_bit = bitand(shifted(valid), neg_shifted);
face_labels(valid) = floor(log2(double(lowest_bit))) + 1;

%% --- Build a colormap for labels ---
unique_labels = unique(face_labels);
num_labels = max(unique_labels) + 1;

% Use HSV for good visual separation; label 0 = dark gray
cmap = uint8(hsv(max(num_labels, 2)) * 255);
cmap = [64 64 64; cmap(2:end, :)];  % label 0 -> gray
if size(cmap,1) < num_labels
    cmap = [cmap; repmat(uint8([128 128 128]), num_labels - size(cmap,1), 1)];
end

face_colors = cmap(face_labels + 1, :);  % Nx3 uint8 RGB

%% --- Write surface PLY with per-face colors ---
fprintf('Saving colored surface mesh to %s ...\n', filename);

fid = fopen(filename, 'w');
nv = size(node, 1);
nf = size(face, 1);

% ASCII header
fprintf(fid, 'ply\n');
fprintf(fid, 'format binary_little_endian 1.0\n');
fprintf(fid, 'element vertex %d\n', nv);
fprintf(fid, 'property float x\n');
fprintf(fid, 'property float y\n');
fprintf(fid, 'property float z\n');
fprintf(fid, 'element face %d\n', nf);
fprintf(fid, 'property list uchar int vertex_indices\n');
fprintf(fid, 'property uchar red\n');
fprintf(fid, 'property uchar green\n');
fprintf(fid, 'property uchar blue\n');
fprintf(fid, 'end_header\n');

% Binary vertices
fwrite(fid, node(:,1:3)', 'float32');

% Binary faces: [count(3), v0, v1, v2, R, G, B] per face
face_idx = int32(face(:, [1 3 2])' - 1);  % 0-based, winding order swap
idx_bytes = typecast(face_idx(:), 'uint8');
idx_bytes = reshape(idx_bytes, 12, nf);

% 1 + 12 + 3 = 16 bytes per face
face_bin = zeros(1, nf * 16, 'uint8');
face_bin(1:16:end) = 3;
for b = 1:12
    face_bin(1+b:16:end) = idx_bytes(b, :);
end
face_bin(14:16:end) = face_colors(:,1);
face_bin(15:16:end) = face_colors(:,2);
face_bin(16:16:end) = face_colors(:,3);

fwrite(fid, face_bin, 'uint8');
fclose(fid);
fprintf('Done: %s (%d vertices, %d faces)\n', filename, nv, nf);

%% --- Optionally save VTU for Paraview (tets + labels) ---
if nargin >= 4 && ~isempty(elem)
    vtu_file = [filename(1:end-4) '.vtu'];
    save_labeled_vtu(vtu_file, node, elem);
end

end


function save_labeled_vtu(filename, node, elem)
% SAVE_LABELED_VTU  Write tetrahedral mesh with region labels as VTK XML (.vtu)
%   Uses binary+base64 encoding for speed and compact file size.
%   Compatible with ParaView, which fully supports tet visualization,
%   slicing, and coloring by region label.

fprintf('Saving labeled tet mesh to %s ...\n', filename);

nv = size(node, 1);
ne = size(elem, 1);

% Decode bitmask labels to integer labels
% bitand(x, -x) isolates lowest set bit, but MATLAB's uint64 negation
% wraps incorrectly, so we use bitcmp (bitwise NOT) to compute two's
% complement: -x == bitcmp(x) + 1
labels = uint64(elem(:, 5));
shifted = bitshift(labels, -1);  % drop bit 0 (background)
decoded = zeros(ne, 1, 'int32');
valid = shifted > 0;
neg_shifted = bitcmp(shifted(valid)) + uint64(1);  % two's complement
lowest_bit = bitand(shifted(valid), neg_shifted);
decoded(valid) = int32(floor(log2(double(lowest_bit)))) + 1;

% Prepare binary arrays
points = single(node(:, 1:3))';      % 3 x nv, column-major -> interleaved xyz
conn = int32(elem(:, 1:4)' - 1);     % 4 x ne, 0-based
offsets = int32((1:ne)' * 4);         % cumulative vertex count per cell
types = uint8(ones(ne, 1) * 10);      % VTK_TETRA = 10

fid = fopen(filename, 'w');

% VTK XML header
fprintf(fid, '<?xml version="1.0"?>\n');
fprintf(fid, '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
fprintf(fid, '  <UnstructuredGrid>\n');
fprintf(fid, '    <Piece NumberOfPoints="%d" NumberOfCells="%d">\n', nv, ne);

% Points
fprintf(fid, '      <Points>\n');
write_b64_array(fid, 'Float32', 3, points);
fprintf(fid, '      </Points>\n');

% Cells
fprintf(fid, '      <Cells>\n');
write_b64_array(fid, 'Int32', 0, conn, 'connectivity');
write_b64_array(fid, 'Int32', 0, offsets, 'offsets');
write_b64_array(fid, 'UInt8', 0, types, 'types');
fprintf(fid, '      </Cells>\n');

% Cell data
fprintf(fid, '      <CellData Scalars="RegionLabel">\n');
write_b64_array(fid, 'Int32', 0, decoded, 'RegionLabel');
fprintf(fid, '      </CellData>\n');

fprintf(fid, '    </Piece>\n');
fprintf(fid, '  </UnstructuredGrid>\n');
fprintf(fid, '</VTKFile>\n');
fclose(fid);

fprintf('Done: %s (%d nodes, %d tets)\n', filename, nv, ne);
end

function write_b64_array(fid, dtype, ncomp, data, name)
% Write a VTK DataArray element with inline base64-encoded binary data.
% VTK inline binary format: 4-byte header (uint32 payload size) + payload,
% all base64-encoded together.
    if nargin < 5, name = ''; end
    
    % Build opening tag
    tag = sprintf('        <DataArray type="%s"', dtype);
    if ~isempty(name), tag = [tag sprintf(' Name="%s"', name)]; end
    if ncomp > 0, tag = [tag sprintf(' NumberOfComponents="%d"', ncomp)]; end
    tag = [tag ' format="binary">\n'];
    fprintf(fid, tag);
    
    % Serialize to raw bytes
    raw = typecast(data(:), 'uint8');
    
    % Prepend uint32 byte count header (VTK requirement for inline binary)
    header = typecast(uint32(numel(raw)), 'uint8');
    payload = [header(:); raw(:)];
    
    % Base64 encode and write
    fprintf(fid, '          %s\n', matlab.net.base64encode(payload));
    fprintf(fid, '        </DataArray>\n');
end
