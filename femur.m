%% Femur Isolation & Morphometrics from NRRD (CT + Slicer Labelmap)
% Requires: Image Processing Toolbox
% Inputs (based on your saved files):
%   1. 'S192803_CT_Scan cropped.nrrd'          - the original CT volume
%   2. 'Segmentation-Segment_1-label.nrrd'      - your femur binary labelmap
%
% NOTE: MATLAB does not have a universal built-in NRRD reader across all
% versions, so a local nrrdread() function is included at the bottom of
% this file - no extra download needed.

clear; clc; close all;

%% ---- 1. SET YOUR PATHS HERE ----
ctPath   = 'E:/bone segmentation/segmentation results/S192803_CT_Scan cropped_1.nrrd';
maskPath = 'E:/bone segmentation/segmentation results/Segmentation_1.nrrd';

%% ---- 2. LOAD CT VOLUME (NRRD) ----
[ctVolume, ctMeta] = nrrdread(ctPath);
ctVolume = double(ctVolume);

% --- HU CALIBRATION FIX ---
% Sanity checks (soft-tissue voxels reading ~1066-1117 instead of the
% expected ~0-100) showed this data is missing the standard DICOM
% RescaleIntercept correction, commonly -1024. Applying it here shifts
% raw scanner values onto the real Hounsfield scale (air=-1000, water=0).
% NOTE: this is an INFERRED correction from internal consistency checks,
% not one confirmed from original scanner metadata - state this plainly
% as a methodological assumption in your report.
HU_OFFSET = -1024;
ctVolume = ctVolume + HU_OFFSET;
fprintf('Applied HU calibration offset: %d\n', HU_OFFSET);

fprintf('CT volume size: %d x %d x %d\n', size(ctVolume,1), size(ctVolume,2), size(ctVolume,3));

% Voxel spacing from NRRD header (mm per voxel, along each axis)
% ctMeta.spacedirections gives the spacing matrix; diagonal = per-axis spacing
spacing = getNrrdSpacing(ctMeta);   % [sx sy sz] in mm
fprintf('Voxel spacing (mm): %.3f x %.3f x %.3f\n', spacing(1), spacing(2), spacing(3));

%% ---- 3. LOAD FEMUR LABELMAP (NRRD) ----
[maskRaw, maskMeta] = nrrdread(maskPath); %#ok<ASGLU>

% --- DIAGNOSTIC: check what's actually in the raw mask data before thresholding ---
fprintf('\n--- Mask Diagnostic ---\n');
fprintf('Mask data class: %s\n', class(maskRaw));
fprintf('Mask unique values (first 10 shown): %s\n', mat2str(unique(maskRaw(1:min(numel(maskRaw),1e6))')));
fprintf('Mask min/max: %g / %g\n', double(min(maskRaw(:))), double(max(maskRaw(:))));
fprintf('Nonzero voxel count: %d\n', nnz(maskRaw));

maskVolume = logical(maskRaw > 0);   % binarize: femur = true, else false

fprintf('Mask volume size: %d x %d x %d\n', size(maskVolume,1), size(maskVolume,2), size(maskVolume,3));

if ~isequal(size(ctVolume), size(maskVolume))
    warning('CT and mask volumes are different sizes! They must match voxel-for-voxel for masking to work correctly.');
end

% --- GEOMETRIC ALIGNMENT CHECK (not just array size) ---
% Same dimensions doesn't guarantee same physical space - origin, spacing,
% and orientation must all match too, or the mask could be silently
% misaligned with the CT despite indexing "working" without error.
fprintf('\n--- Geometry Alignment Check ---\n');
ctSpacing   = getNrrdSpacing(ctMeta);
maskSpacing = getNrrdSpacing(maskMeta);
fprintf('CT spacing:   %s mm\n', mat2str(round(ctSpacing,4)));
fprintf('Mask spacing: %s mm\n', mat2str(round(maskSpacing,4)));
if isfield(ctMeta,'space_origin'), fprintf('CT origin:   %s\n', ctMeta.space_origin); end
if isfield(maskMeta,'space_origin'), fprintf('Mask origin: %s\n', maskMeta.space_origin); end
if max(abs(ctSpacing - maskSpacing)) > 1e-3
    warning('CT and mask spacing differ by more than 0.001mm - they may not be co-registered.');
end
if isfield(ctMeta,'space_origin') && isfield(maskMeta,'space_origin')
    ctOrigin   = sscanf(ctMeta.space_origin, '(%f,%f,%f)');
    maskOrigin = sscanf(maskMeta.space_origin, '(%f,%f,%f)');
    originDiff = max(abs(ctOrigin - maskOrigin));
    fprintf('Maximum origin difference: %.12g mm\n', originDiff);
    if originDiff > 1e-3
        warning('CT and mask origins differ by >0.001mm - verify these volumes are truly aligned, not just same-sized.');
    else
        fprintf('Origins match within tolerance - CT and mask are properly co-registered.\n');
    end
end

% --- HU CALIBRATION NOTE ---
% NRRD does not carry DICOM RescaleSlope/RescaleIntercept tags. If this
% file was exported from Slicer after loading DICOM, Slicer typically
% already bakes the rescale into the voxel values (so raw NRRD values
% ARE true HU). But this is an assumption, not a guarantee - if your
% source wasn't DICOM, or a custom import path was used, verify against
% a known landmark (e.g. air outside the body should read ~-1000 HU,
% water/soft tissue ~0-100 HU) before treating these as calibrated HU.
if min(ctVolume(:)) >= 0
    fprintf('\nNOTE: CT min value is %.0f (not negative). True calibrated HU air is ~-1000.\n', min(ctVolume(:)));
    fprintf('This suggests either: (a) your crop excluded all air/background, or (b) values are offset/uncalibrated.\n');
    fprintf('Worth confirming before treating these as standard HU.\n');
end

if nnz(maskVolume) == 0
    error(['Mask contains ZERO foreground voxels after thresholding (maskVolume > 0). ' ...
           'Check the diagnostic printout above - if min/max are both 0, the label file ' ...
           'itself has no segmented data, or is being read incorrectly. ' ...
           'Verify the file in Slicer/a NRRD viewer first.']);
end

%% ---- 4. MASK OUT EVERYTHING EXCEPT THE FEMUR ----
femurOnly = ctVolume .* maskVolume;
femurOnly(~maskVolume) = NaN;   % background -> NaN so stats ignore it

%% ---- 5. CROP TO BOUNDING BOX ----
stats3D = regionprops3(maskVolume, 'BoundingBox');
bb = round(stats3D.BoundingBox(1,:));   % regionprops3 order: [y_start x_start z_start y_width x_width z_width]

yRange = bb(2):(bb(2)+bb(5)-1);
xRange = bb(1):(bb(1)+bb(4)-1);
zRange = bb(3):(bb(3)+bb(6)-1);

femurCropped = femurOnly(yRange, xRange, zRange);
maskCropped  = maskVolume(yRange, xRange, zRange);

%% ---- 6. BASIC MORPHOMETRIC STATS ----
props = regionprops3(maskVolume, 'Volume', 'BoundingBox', 'PrincipalAxisLength', 'Centroid');
voxelVol_mm3 = spacing(1) * spacing(2) * spacing(3);
totalVolume_mm3 = props.Volume(1) * voxelVol_mm3;

fprintf('\n--- Femur Morphometry ---\n');
fprintf('Voxel count: %d voxels\n', props.Volume(1));
fprintf('Physical volume: %.1f mm^3 (%.2f cm^3)\n', totalVolume_mm3, totalVolume_mm3/1000);
fprintf('Bounding box (voxels): %s\n', mat2str(round(props.BoundingBox(1,:))));
fprintf('Principal axis lengths (voxels): %s\n', mat2str(round(props.PrincipalAxisLength(1,:))));

%% ---- 7. HOUNSFIELD UNIT (DENSITY) DISTRIBUTION ----
huValues = femurCropped(~isnan(femurCropped));

figure('Name','Femur HU Distribution');
histogram(huValues, 50);
xlabel('Hounsfield Units (HU)');
ylabel('Voxel Count');
title('Bone Density Distribution - Femur Only');
grid on;

fprintf('\n--- HU Statistics (Full Distribution) ---\n');
fprintf('Mean HU: %.1f\n', mean(huValues));
fprintf('Median HU: %.1f\n', median(huValues));
fprintf('Std HU: %.1f\n', std(huValues));
fprintf('Min/Max HU: %.1f / %.1f\n', min(huValues), max(huValues));

p = prctile(huValues, [5 25 75 95]);
fprintf('5th percentile: %.1f\n', p(1));
fprintf('25th percentile (Q1): %.1f\n', p(2));
fprintf('75th percentile (Q3): %.1f\n', p(3));
fprintf('95th percentile: %.1f\n', p(4));
fprintf('Interquartile range (IQR): %.1f\n', p(3)-p(2));

fprintf('\n--- Density Band Breakdown ---\n');
bands = [-inf 300; 300 600; 600 1000; 1000 1500; 1500 2000; 2000 inf];
bandLabels = {'<300 HU','300-600 HU','600-1000 HU','1000-1500 HU','1500-2000 HU','>2000 HU'};
for i = 1:size(bands,1)
    pct = 100 * sum(huValues >= bands(i,1) & huValues < bands(i,2)) / numel(huValues);
    fprintf('%-15s: %5.2f%%\n', bandLabels{i}, pct);
end

% --- SPATIAL LOCALIZATION OF HIGH-HU VOXELS (investigate the >2000 HU peak) ---
% Find where in the volume these high-density voxels actually sit, so we
% can tell if it's real dense cortical bone (should cluster at the outer
% shell/mid-shaft), a segmentation leak (would show at mask edges/outside
% expected bone shape), or scattered noise (would be randomly distributed).
highHU_mask = femurCropped > 2000;
if nnz(highHU_mask) > 0
    [ri, ci, zi] = ind2sub(size(femurCropped), find(highHU_mask));
    fprintf('\n--- >2000 HU Voxel Locations ---\n');
    fprintf('Count: %d voxels (%.2f%% of femur)\n', numel(zi), 100*numel(zi)/numel(huValues));
    fprintf('Z-slice range (within cropped volume): %d to %d (out of %d total slices)\n', min(zi), max(zi), size(femurCropped,3));
    fprintf('If this range clusters narrowly (e.g. only a few slices), investigate those slices specifically for artifacts.\n');
    fprintf('If spread evenly across most slices, likely genuine dense cortical bone (expected at the outer shell).\n');

    % Quick visual: show where these voxels sit relative to the whole cropped mask
    figure('Name','High-HU Voxel Locations (>2000 HU)');
    highHUCountPerSlice = squeeze(sum(sum(highHU_mask,1),2));
    totalCountPerSlice = squeeze(sum(sum(maskCropped,1),2));
    plot(1:numel(highHUCountPerSlice), highHUCountPerSlice, 'r-', 'LineWidth', 1.5); hold on;
    plot(1:numel(totalCountPerSlice), totalCountPerSlice/20, 'k--'); % scaled down for visibility on same axis
    legend('>2000 HU voxel count per slice', 'Total femur voxel count per slice (/20, for scale)');
    xlabel('Slice index (within cropped volume)');
    ylabel('Voxel count');
    title('Where the >2000 HU Voxels Are Located Along the Bone');
    grid on;
else
    fprintf('\nNo voxels above 2000 HU found in this run.\n');
end

%% ---- 8. CROSS-SECTIONAL AREA vs POSITION ALONG SHAFT ----
numSlices = size(maskCropped, 3);
crossSectionalArea = zeros(numSlices,1);
pixelArea_mm2 = spacing(1) * spacing(2);

for z = 1:numSlices
    crossSectionalArea(z) = sum(sum(maskCropped(:,:,z))) * pixelArea_mm2;
end

sliceThickness = spacing(3);
positionAlongShaft_mm = (0:numSlices-1) * sliceThickness;

figure('Name','Cross-Sectional Area Along Shaft');
plot(positionAlongShaft_mm, crossSectionalArea, 'LineWidth', 1.5);
xlabel('Position Along Shaft (mm)');
ylabel('Cross-Sectional Area (mm^2)');
title('Femur Cross-Sectional Area Profile');
grid on;

%% ---- 8b. ANATOMICAL REGIONAL BREAKDOWN ----
% Orientation established from NRRD header: 'space directions' has a
% POSITIVE Z-spacing (0,0,0.6) under 'left-posterior-superior' (LPS)
% convention, meaning increasing slice index = increasing SUPERIOR
% direction = toward the femoral head. So:
%   low slice index  (position ~0mm)   = condyles (distal/inferior)
%   high slice index (position ~max mm) = femoral head/neck (proximal/superior)
%
% If your bone renders upside-down relative to this in Slicer, set
% FLIP_ORIENTATION = true below to reverse the region assignment.
FLIP_ORIENTATION = false;
 
% Region boundaries as fraction of total bone length (adjust if you want
% to match specific anatomical landmarks visible in your cross-sectional
% area plot instead of even percentage splits)
regionBoundaryFrac = [0 0.10 0.25 0.75 0.90 1.0];
regionNames = {'Condyles (distal)', 'Distal metaphysis', 'Midshaft', ...
                'Proximal metaphysis', 'Femoral head/neck (proximal)'};

if FLIP_ORIENTATION
    regionNames = fliplr(regionNames);
end

sliceBoundaries = round(regionBoundaryFrac * numSlices);
sliceBoundaries(1) = 1;
sliceBoundaries(end) = numSlices;

CORTICAL_THRESHOLD = 600; % HU cutoff for "cortical/high-density" fraction, matching earlier density-band boundary

fprintf('\n\n========== ANATOMICAL REGIONAL BREAKDOWN ==========\n');
fprintf('%-30s %8s %8s %8s %8s %8s %8s %10s %10s %10s\n', ...
    'Region', 'Mean', 'Median', 'P5', 'P25', 'P75', 'P95', 'Cortical%', 'Vol(cm3)', 'MeanCSA');

regionResults = struct();
for i = 1:5
    sliceStart = sliceBoundaries(i) + (i>1);
    sliceEnd = sliceBoundaries(i+1);
    if sliceStart > sliceEnd, sliceStart = sliceEnd; end
    sliceRange = sliceStart:sliceEnd;

    huRegion = femurCropped(:,:,sliceRange);
    huRegion = huRegion(~isnan(huRegion));

    maskRegion = maskCropped(:,:,sliceRange);
    volRegion_mm3 = sum(maskRegion(:)) * voxelVol_mm3;

    meanCSA = mean(crossSectionalArea(sliceRange));
    corticalPct = 100 * sum(huRegion > CORTICAL_THRESHOLD) / numel(huRegion);
    pReg = prctile(huRegion, [5 25 75 95]);

    fprintf('%-30s %8.1f %8.1f %8.1f %8.1f %8.1f %8.1f %9.2f%% %10.2f %10.1f\n', ...
        regionNames{i}, mean(huRegion), median(huRegion), pReg(1), pReg(2), pReg(3), pReg(4), ...
        corticalPct, volRegion_mm3/1000, meanCSA);

    regionResults.(matlab.lang.makeValidName(regionNames{i})) = struct( ...
        'sliceRange', sliceRange, 'huValues', huRegion, 'mean', mean(huRegion), ...
        'median', median(huRegion), 'percentiles', pReg, 'corticalPct', corticalPct, ...
        'volume_mm3', volRegion_mm3, 'meanCSA_mm2', meanCSA);
end
fprintf('=====================================================\n');

% Bar chart comparing mean HU across regions
figure('Name','Regional HU Comparison');
regionMeans = cellfun(@(n) regionResults.(matlab.lang.makeValidName(n)).mean, regionNames);
bar(categorical(regionNames, regionNames), regionMeans);
ylabel('Mean HU');
title('Mean HU by Anatomical Region');
grid on;

%% ---- 9. 3D RENDER OF THE SEGMENTED FEMUR (isosurface) ----
figure('Name','3D Femur Reconstruction');
fv = isosurface(maskCropped, 0.5);   % 0.5 = threshold between 0 (background) and 1 (femur)

% Scale the mesh vertices to real-world mm using voxel spacing
% (isosurface works in voxel index space by default, so we rescale)
fv.vertices = fv.vertices .* spacing([2 1 3]); % note x/y swap due to MATLAB array convention

p = patch(fv);
set(p, 'FaceColor', [0.85 0.75 0.6], 'EdgeColor', 'none');
daspect([1 1 1]);
view(3);
axis tight;
camlight('headlight');
lighting gouraud;
material dull;
xlabel('X (mm)'); ylabel('Y (mm)'); zlabel('Z (mm)');
title('3D Reconstructed Femur (from MATLAB segmentation mask)');
rotate3d on;

%% ---- 10. SAVE RESULTS ----
save('femur_morphometry_results.mat', 'huValues', 'crossSectionalArea', 'props', 'totalVolume_mm3', 'regionResults');
fprintf('\nResults saved to femur_morphometry_results.mat\n');


%% ================= LOCAL FUNCTIONS =================

function [data, meta] = nrrdread(filename)
% Minimal NRRD reader (detached or attached header, raw/gzip encoding)
    fid = fopen(filename, 'rb');
    if fid < 0, error('Could not open file %s', filename); end

    meta = struct();
    line = fgetl(fid);
    if ~strncmp(line, 'NRRD', 4)
        error('Not a valid NRRD file.');
    end

    % Read header key/value lines until blank line
    while true
        line = fgetl(fid);
        if isempty(line), break; end
        if line(1) == '#', continue; end
        colonIdx = strfind(line, ':');
        if isempty(colonIdx), continue; end
        key = strtrim(line(1:colonIdx(1)-1));
        val = strtrim(line(colonIdx(1)+1:end));
        key = strrep(key, ' ', '_');
        key = strrep(key, '-', '');
        meta.(matlab.lang.makeValidName(key)) = val;
    end

    % Parse dimensions and type
    sizes = sscanf(meta.sizes, '%d')';
    typeMap = struct( ...
        'signedchar','int8', 'int8','int8', 'int8_t','int8', ...
        'uchar','uint8', 'unsignedchar','uint8', 'uint8_t','uint8', ...
        'short','int16', 'shortint','int16', 'signedshort','int16', 'signedshortint','int16', 'int16_t','int16', ...
        'ushort','uint16', 'unsignedshort','uint16', 'unsignedshortint','uint16', 'uint16_t','uint16', ...
        'int','int32', 'signedint','int32', 'int32_t','int32', ...
        'uint','uint32', 'unsignedint','uint32', 'uint32_t','uint32', ...
        'longlong','int64', 'signedlonglong','int64', 'longlongint','int64', 'int64_t','int64', ...
        'ulonglong','uint64', 'unsignedlonglong','uint64', 'unsignedlonglongint','uint64', 'uint64_t','uint64', ...
        'float','single', 'double','double');

    % Normalize: lowercase, strip all whitespace and hyphens so
    % "unsigned short" / "unsigned-short" / "unsignedshort" all match
    typeStr = lower(strtrim(meta.type));
    typeStr = regexprep(typeStr, '[\s\-]', '');

    if isfield(typeMap, typeStr)
        matType = typeMap.(typeStr);
    else
        error('Unrecognized NRRD data type: "%s". Please tell me this value so I can add it to the type map.', meta.type);
    end

    encoding = lower(strtrim(meta.encoding));

    rawData = fread(fid, Inf, 'uint8=>uint8');
    fclose(fid);

    if strcmp(encoding, 'gzip') || strcmp(encoding, 'gz')
        rawData = gunzipBytes(rawData);
    elseif ~strcmp(encoding, 'raw')
        error('Unsupported NRRD encoding: %s', encoding);
    end

    data = typecast(rawData, matType);
    data = reshape(data, sizes);
end

function out = gunzipBytes(bytes)
% Decompress gzip byte stream using MATLAB's built-in gunzip (file-based).
% More reliable for large files than manual Java streaming, which can
% silently truncate/fail on large buffers depending on JVM heap settings.
    tmpGzFile = [tempname, '.gz'];
    fid = fopen(tmpGzFile, 'wb');
    fwrite(fid, bytes, 'uint8');
    fclose(fid);

    outFiles = gunzip(tmpGzFile, tempdir);

    fid2 = fopen(outFiles{1}, 'rb');
    out = fread(fid2, Inf, 'uint8=>uint8');
    fclose(fid2);

    % Clean up temp files
    delete(tmpGzFile);
    delete(outFiles{1});
end

function spacing = getNrrdSpacing(meta)
% Extract per-axis voxel spacing (mm) from NRRD header metadata
% Field names get 'space'->'_' substitution applied during header parsing,
% so check both underscored and concatenated variants to be safe.
    if isfield(meta, 'spacings')
        vals = sscanf(meta.spacings, '%f');
        spacing = vals(:)';
    elseif isfield(meta, 'space_directions')
        spacing = spacingFromDirections(meta.space_directions);
    elseif isfield(meta, 'spacedirections')
        spacing = spacingFromDirections(meta.spacedirections);
    else
        warning('Could not find spacing info in NRRD header - defaulting to 1mm isotropic. Available header fields: %s', strjoin(fieldnames(meta), ', '));
        spacing = [1 1 1];
    end
end

function spacing = spacingFromDirections(str)
% Format like: (1,0,0) (0,1,0) (0,0,1) or "none (1,0,0) (0,1,0) (0,0,1)"
% (leading "none" appears when an extra non-spatial axis is present)
    nums = regexp(str, '[-\d\.]+', 'match');
    nums = str2double(nums);
    nums = nums(~isnan(nums));   % drop anything that wasn't actually numeric
    nums = nums(1:9);            % keep only the 3x3 spatial part if extra values present
    M = reshape(nums, 3, 3)';
    spacing = sqrt(sum(M.^2, 2))';
end