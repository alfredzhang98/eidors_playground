function out = mk_thorax_model_bp3d(varargin)
%MK_THORAX_MODEL_BP3D Torso model based on BodyParts3D
%
% FMDL = MK_THORAX_MODEL_BP3D() 
% Returns a torso fwd_model with no electrodes or organs
%
% FMDL = MK_THORAX_MODEL_BP3D(elec_pos, elec_shape)
% FMDL = MK_THORAX_MODEL_BP3D(elec_pos, elec_shape, maxh)
% Allows specifying electrode configuration, see PLACE_ELEC_ON_SURF
%
% FMDL = MK_THORAX_MODEL_BP3D(stimpat)
% A predefined model with round electrodes (radius 1 cm). Accepts a string 
% parameter to define a stim pattern:
%    '2x16_planar'   - 2 rings of 16 electrodes, planar stimulation
%    '2x16_odd-even' - 2 rings, each stimulation is cross-plane 
%    '2x16_square'   - 2 rings, every second stimulation is cross-plane
%    '2x16_adjacent' - like square, but with adjacent stimulation
%    '1x32_ring'     - single ring of 32 electrodes
%
% IMG = MK_THORAX_MODEL_BP3D(... , organspec)
% Returns a torso model containing organ shapes as an image structure. Can
% be used with all the above electrode defintions. 
% ORGANSPEC must be one of:
%    'fine'          - a relatively dense model, ~7 million elems
%    'coarse'        - a coarser model, *NOT IMPLEMENTED YET*
% 
%
% The 'fine' model contains the following materials:
%    #   mat_names         Description
%   -----------------------------------------------------------------------
%     1	Body             	Volume between the skin and the internal organs
%     2	LungLeft         	Left lung
%     3	LungRight        	Right lung
%     4	HeartWall        	Heart wall
%     5	Aorta            	Aorta (wall)
%     6	SVC              	Superior Vena Cava (wall)
%     7	IVC              	Inferior Vena Cava (wall)
%     8	PA               	Pulmonary Artery (wall)
%     9	LSPV             	Left Superior Pulmonary Vein (wall)
%    10	LIPV             	Left Inferior Pulmonary Vein (wall)
%    11	RSPV             	Right Superior Pulmonary Vein (wall)
%    12	RIPV             	Right Inferior Pulmonary Vein (wall)
%    13	BloodHeartLeft   	Blood in left atrium and left ventricle
%    14	BloodHeartRight  	Blood in right atrium and right ventricle
%    15	BloodAorta       	Aorta (lumen)
%    16	BloodSVC         	Superior Vena Cava (lumen)
%    17	BloodIVC         	Inferior Vena Cava (lumen)
%    18	BloodPA          	Pulmonary Artery (lumen)
%    19	BloodLSPV        	Left Superior Pulmonary Vein (lumen)
%    20	BloodLIPV        	Left Inferior Pulmonary Vein (lumen)
%    21	BloodRSPV        	Right Superior Pulmonary Vein (lumen)
%    22	BloodRIPV        	Right Inferior Pulmonary Vein (lumen)
%
%
% CITATION_REQUEST:
% AUTHOR: Grychtol B
% TITLE: Torso Organ Meshes
% JOURNAL: Zenodo 
% YEAR: 2024
% DOI: 10.5281/zenodo.11047831
%
% LICENSE INFORMATION:
% The model generated by this function is a derivative of the <a href=
% "matlab: web https://https://doi.org/10.5281/zenodo.11047830 -browser"
% >"Torso Organ Meshes"</a> 
% by Bartłomiej Grychtol, itself a derivative of the <a href="matlab: web 
% https://dbarchive.biosciencedbc.jp/en/bodyparts3d -browser">BodyParts3D</a> database by 
% <a href="matlab: web https://biosciencedbc.jp/ -browser"
% >The Database Center for Life Science</a>, used under the Creative Commons 
% Attribution-ShareAlike 2.1 Japan (<a href="matlab: web 
% http://creativecommons.org/licenses/by-sa/2.1/jp/deed.en_US 
% -browser">CC BY-SA 2.1 JP</a>) license.
% The model must be licensed under CC BY-SA 2.1 JP and attributed
% correctly.
%
% See also MK_THORAX_MODEL, PLACE_ELEC_ON_SURF
%
% For usage examples: pubs/conferences/EIT2024/torso-mesh-bp3d/mk_figures.m

% (C) 2024 Bartek Grychtol. License: GPL v2 or v3
% $Id: mk_thorax_model_bp3d.m 7137 2024-12-29 18:42:23Z aadler $

if nargin>=1 && strcmp(varargin{1},'UNIT_TEST'); do_unit_test; return; end

citeme(mfilename)

if nargin == 0 % get a torso with no electrodes
   out = get_torso; 
end
organspec = '';
copyright = '';
if nargin > 1 && ischar(varargin{end})
   organspec = varargin{end};
   varargin(end) = [];
end

switch length(varargin)
   case 0
      out = mk_organ_image(out, organspec); % for the sake of completeness 
   case 1
      if isempty(varargin{1}) 
         out = get_torso;
         if ~isempty(organspec)
            out = mk_organ_image(out, organspec);
         end
      elseif ischar(varargin{1}) 
         copyright = '(C) EIDORS Project';
         out = build_predef_elecs(varargin{1}, organspec);
      else
         error('EIDORS:WrongInput','Electrode specification not understood');
      end
   case {2,3}
      out = mk_thorax_model('bp3d',varargin{:});
      out = mk_organ_image(out, organspec);
   otherwise
      error('EIDORS:WrongInput', 'Too many inputs');
end

out = eidors_obj('set', out); % sort fields

if ~isempty(copyright)
   out.copyright = copyright;
end

out.attribution = ...
   ['Modified from "Torso Organ Meshes" ' ...
    '(https://doi.org/10.5281/zenodo.11047831) by Bartłomiej Grychtol ' ...
    '/ A derivative of "BodyParts3D" '...
    '(https://dbarchive.biosciencedbc.jp/en/bodyparts3d) by '...
    'The Database Center for Life Science.'];

out.license   =  ['Creative Commons Attribution-ShareAlike 2.1 Japan '...;
       '(http://creativecommons.org/licenses/by-sa/2.1/jp/deed.en_US)'];

if strcmp(out.type, 'image')
   try out.fwd_model.copyright = out.copyright; end
   out.fwd_model.attribution = out.attribution;
   out.fwd_model.license = out.license;
end

ws = warning('off', 'backtrace');
warning('EIDORS:License', ...
   ['The bp3d torso model is licensed under CC BY-SA 2.1 JP.\n' ...
    'Please make sure to attribute correctly. ' ...
    'See the attribution field for details.']);
warning(ws.state, 'backtrace');


function out = get_torso()
contrib = 'bp3d-torso'; file =  'bp3d_torso.mat';
load(get_contrib_data(contrib,file),'torso');
eidors_msg('Building bp3d body surface', 2);
out = gmsh_stl2tet(torso);
out = fix_boundary(out);
out.name = 'BP3D Torso Model';


function out = build_predef_elecs(stimpat, organspec)
cache_path = get_cache_path;
fname = 'bp3d_predef_torso';
suffix = '';
if ~isempty(organspec)
   fname = [fname '_' organspec];
   suffix = [' - ' organspec];
end
name = ['BP3D Torso Model' suffix];
fpath = [cache_path filesep fname '.mat'];
if exist(fpath,'file')
   eidors_msg('@@ Using stored model.', 3);
   load(fpath, 'out');
else
   eth16 = 360*cumsum([0.2 0.4*ones(1,7) 0.5 0.4*ones(1,7)])/6.5 - 90; eth16 = eth16';
   eth32 = 360*cumsum([0.1 0.2*ones(1,15) 0.25 0.2*ones(1,15)])/6.45 - 90; eth32 = eth32';
   ep = eth16; ep(:,2) = 1175;
   ep(17:48,1) = eth32; ep(17:48,2) = 1212.5;
   ep(49:64,1) = eth16; ep(49:64,2) = 1250;
   out = mk_thorax_model('bp3d',ep,[10 0 1],20);
   out.name = name;
   out = mk_organ_image(out, organspec);
   save(fpath,'out');
end
out = set_predef_stim(out, stimpat);
out.name = [out.name ' - ' stimpat];


function out = mk_organ_image(out, organspec)
if isempty(organspec), return, end
organs = get_organ_models(organspec);
out = eidors_cache(@add_organs,{out, organs});
out = set_organ_values(out);
out.name = 'BP3D Torso Image';


function out = add_organs(torso, organs)
tri.vertices = torso.nodes;
tri.faces = torso.boundary;
N = size(tri.vertices, 1);
tri.vertices = [tri.vertices; organs.nodes];
tri.faces = [tri.faces; N + organs.boundary(:,[1 3 2])];

eidors_msg('Re-meshing body volume', 2);
shell = gmsh_stl2tet(tri,[],[],4); % lots of optimization
print_mdl_info(shell, 'Body');

eidors_msg('Restoring electrodes', 2);
if isfield(torso, 'electrode') && ~isempty(torso.electrode)
   for i = 1:numel(torso.electrode)
      electrode(i).nodes = torso.nodes(torso.electrode(i).nodes,:);
   end
   shell = restore_electrodes(shell, electrode, 0.01);
end
organs.boundary = organs.show_boundary;
organs = rmfield(organs, 'show_boundary');

eidors_msg('Merging body and organ volumes', 2);
out = merge_meshes(shell, organs);
out.mat_names = ['Body'; organs.mat_names];
 
eidors_msg('Removing trachea', 2);
idx = strcmp(out.mat_names,'Trachea');
out = remove_elems(out, out.mat_idx{idx});

% nicely re-order materials
idx = [1 5 10 2  17:2:21 15 6 8 11 13 3 4 18:2:22 16 7 9 12 14];
out.mat_idx = out.mat_idx(idx);
out.mat_names = out.mat_names(idx);

out = eidors_obj('set', out); % order fields

print_mdl_info(out, 'Final model');

function img = set_organ_values(torso)
img = mk_image(torso, 1, 'bp3d-torso');
for i = 2:numel(torso.mat_idx)
    if strcmp(torso.mat_names{i},'HeartWall')
        img.elem_data(torso.mat_idx{i}) = 1.5;
%     elseif strcmp(torso.mat_names{i},'Trachea')
%         img.elem_data(torso.mat_idx{i}) = NaN;
    elseif ~isempty(strfind(torso.mat_names{i}, 'Blood'))
        img.elem_data(torso.mat_idx{i}) = 2;
    elseif ~isempty(strfind(torso.mat_names{i}, 'Lung'))
        img.elem_data(torso.mat_idx{i}) = .5;
    else % vessel
        img.elem_data(torso.mat_idx{i}) = .8;
    end
end
img.calc_colours.ref_level = 1;

function mm = restore_electrodes(mm, electrode, z_contact, thresh)
if nargin < 4
   thresh = eps;
end
surf_node_idx = unique(find_boundary(mm));
surf_nodes = mm.nodes(surf_node_idx,:);
% ignore nodes close to the bottom surface
rng = max(surf_nodes(:,3)) - min(surf_nodes(:,3));
idx = surf_nodes(:,3) > min(surf_nodes(:,3)) + 0.02 * rng;
surf_node_idx = surf_node_idx(idx);
surf_nodes = surf_nodes(idx,:);
for e = 1:numel(electrode)
   % if we got this far, we should have enough memory for vectorized calcs
   dist = (surf_nodes(:,1) - electrode(e).nodes(:,1)').^2;
   for d = 2:3
      dist = dist + (surf_nodes(:,d) - electrode(e).nodes(:,d)').^2;
   end
%    dist = sqrt(dist); % takes a long time and changes nothing
   [val, pos] = min(dist);
   if any(val > thresh)
      warning('Some nodes of electrode %d seem too far. Furthest distance: %f',e,sqrt(max(val)));
   end
   mm.electrode(e).nodes = surf_node_idx(pos);
   mm.electrode(e).z_contact = z_contact;
end


function organs = get_organ_models(organspec)
switch organspec
   case 'fine'
      fhandle = @mk_fine_organs;
   case 'coarse'
      error('EIDORS:NotImplemented', ['Coarse organ model is not '...
            'implemented yet.'])
   otherwise
      error('EIDORS:WrongInput', 'Organ specification not understood.');
end
cache_path = get_cache_path;
fname = sprintf('bp3d_organs_%s.mat', organspec);
fpath = [cache_path filesep fname];
if ~exist(fpath, 'file')
   eidors_msg('This takes a REALLY long time');
   organs = feval(fhandle);
   save(fpath, 'organs');
else
   load(fpath, 'organs');
end


function mdl = mk_fine_organs

start_time = now(); 
organs = mk_individual_organs;
mdl = assemble_organs(organs);
mdl.show_boundary = mdl.boundary;
eidors_msg('Recalculating boundary',2);
mdl = fix_boundary(mdl);

%  cashed results are not useful outside
eidors_cache('clear_old', start_time);


function mdl = assemble_organs(organs)
pcs = struct();

eidors_msg('Merging vessels', 2);
vessels = {'Aorta','IVC','SVC','PA', 'LSPV', 'LIPV','RSPV','RIPV'};
for i = 1:numel(vessels)
    v = vessels{i};
    disp(v)
    pcs.(v) = merge_meshes(organs.(v), organs.(['Blood' v]),1e-2);
    pcs.(v).mat_names = {v; ['Blood' v]};
end

eidors_msg('Merging left side', 2);
left = merge_meshes(organs.LungLeft, ...            
                    pcs.LSPV, ...
                    pcs.LIPV, ...
                    1e-2);
left.mat_names = vertcat({'LungLeft'},pcs.LSPV.mat_names, pcs.LIPV.mat_names);

eidors_msg('Merging right side', 2);
right = merge_meshes(organs.LungRight, ...            
                     pcs.RSPV, ...
                     pcs.RIPV, ...
                     2e-1); % need a larger threshold due to mesh mismatch
right.mat_names = vertcat({'LungRight'},pcs.RSPV.mat_names, pcs.RIPV.mat_names);

eidors_msg('Merging heart', 2);
heart = merge_meshes(organs.HeartWall, ...
                     organs.BloodLeftHeart, ...
                     organs.BloodRightHeart, ...
                     1e-2);
heart.mat_names = {'HeartWall';'BloodHeartLeft'; 'BloodHeartRight'};

eidors_msg('Merging heart, left and right', 2);
mdl = merge_meshes(heart, left, right, 1e-2);
mdl.mat_names = vertcat(heart.mat_names, left.mat_names, right.mat_names);

eidors_msg('Merging remaining vessels', 2);
mdl = merge_meshes(mdl, ...
                  pcs.PA, ...
                  pcs.Aorta, ...
                  pcs.SVC, ...
                  pcs.IVC, ...
                  1e-2);
mdl.mat_names = vertcat(mdl.mat_names, pcs.PA.mat_names,pcs.Aorta.mat_names, ...
                        pcs.SVC.mat_names, pcs.IVC.mat_names);

eidors_msg('Merging trachea', 2);
mdl = merge_meshes(mdl, organs.Trachea, 5e-2);
mdl.mat_names{end+1} = 'Trachea';


function organs = mk_individual_organs
org_path = get_contrib_data('bp3d-torso','bp3d_organs.mat');
srf = load(org_path);

eidors_msg('Meshing vessels', 2);
vessels = {'Aorta','IVC','SVC','PA', 'LSPV', 'LIPV','RSPV','RIPV'};
for i = 1:numel(vessels)
    v = vessels{i};
    disp(v)
    [vessel, blood] = eidors_cache(@mk_vessel, srf.(v));
    organs.(v) = vessel;
    b = sprintf('Blood%s', v);
    organs.(b) = blood;
    print_mdl_info(vessel, v);
    print_mdl_info(blood, b);
end

eidors_msg('Meshing trachea and lungs', 2);
organs.Trachea = gmsh_stl2tet(srf.Trachea,[],[],1);
print_mdl_info(organs.Trachea, 'Trachea');
organs.LungLeft = gmsh_stl2tet(srf.LungLeft,[],[],2);
print_mdl_info(organs.LungLeft, 'LungLeft');
organs.LungRight = gmsh_stl2tet(srf.LungRight,[],[],2);
print_mdl_info(organs.LungRight, 'LungRight');

eidors_msg('Meshing left heart', 2);
left = eidors_cache(@mk_blood_left_surf,srf);
organs.BloodLeftHeart = gmsh_stl2tet(left, [3, 20],[],4);
print_mdl_info(organs.BloodLeftHeart, 'BloodLeftHeart');

eidors_msg('Meshing right heart', 2);
right = eidors_cache(@mk_blood_right_surf,srf);
organs.BloodRightHeart = gmsh_stl2tet(right, [3, 20],[],4);
print_mdl_info(organs.BloodRightHeart, 'BloodRightHeart');

eidors_msg('Meshing heart wall', 2);
heart = eidors_cache(@mk_heart_surf, srf);
organs.HeartWall = gmsh_stl2tet(heart,[],4);
print_mdl_info(organs.HeartWall, 'HeartWall');


function right = mk_blood_right_surf(srf)
fn = @(x) flip_normals(x);
right = merge_meshes(srf.BloodRight, ...
                     fn(srf.SVC.cap), ...
                     fn(srf.IVC.cap), ...
                     fn(srf.PA.cap));


function left = mk_blood_left_surf(srf)
fn = @(x) flip_normals(x);
left = merge_meshes(srf.BloodLeft, ...
                    fn(srf.Aorta.cap), ...
                    fn(srf.LSPV.cap), ...
                    fn(srf.LIPV.cap), ...
                    fn(srf.RSPV.cap), ...
                    fn(srf.RIPV.cap), ...
                    1e-2);


function [heart, right, left] = mk_heart_surf(srf)
right = merge_meshes(srf.BloodRight, ...
                     srf.SVC.ring, ...
                     srf.IVC.ring, ...
                     srf.PA.ring, ...
                     1e-2);
right = flip_normals(right);

left = merge_meshes(srf.BloodLeft, ...
                    srf.Aorta.ring, ...
                    srf.LSPV.ring, ...
                    srf.LIPV.ring, ...
                    srf.RSPV.ring, ...
                    srf.RIPV.ring, ...
                    1e-2);
left = flip_normals(left);

heart = merge_meshes(bd(left), bd(right), bd(srf.HeartShell));


function [V, B] = mk_vessel(vstruct)
msh = merge_meshes(vstruct.outer, vstruct.ring, ...
                    flip_normals(vstruct.inner),1e-2);
V = gmsh_stl2tet(msh, 1,[],3); % max 1 mm edge (vessels are 2 mm thick
msh = merge_meshes(vstruct.cap, vstruct.inner);
B = gmsh_stl2tet(msh,[],1);


function mdl = flip_normals(mdl)
mdl.elems = mdl.elems(:,[1 3 2]);


function mdl = bd(mdl)
mdl.boundary = mdl.elems;


function print_mdl_info(mdl, name)
fprintf('=> %20s: %6d nodes, %7d elems\n',name, num_nodes(mdl), num_elems(mdl));


function out = get_cache_path
global eidors_objects
out = [eidors_objects.model_cache filesep 'bp3d-torso'];
if ~exist(out, 'dir')
   mkdir(out);
end

function do_unit_test
   torso = mk_thorax_model_bp3d('2x16_odd-even','fine');
   expect = [206.019760131836   45.109825134277   1641.201293945312];
   unit_test_cmp('extent',max(torso.fwd_model.nodes), expect,1e-10);
   unit_test_cmp('#elems', num_elems(torso), 6894847);
   unit_test_cmp('#nodes', num_nodes(torso), 1193751);
