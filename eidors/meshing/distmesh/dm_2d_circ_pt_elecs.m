function fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, spacing);
% DM_2D_CIRC_PT_ELECS: Create circle mesh (or radius 1) refined with electrodes
% fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, spacing);
%     at points on the electrodes
% fmdl = dm_2d_circ_pt_elecs( elec_pts, pfix, ...
%               [base_spacing, refine_ratio, gradient]);
% elec_pts = cell fcn of N x [x,y,{z}] for each electrode
%    normally only two points are specified (at start and end of electrode)
%    eg elec_pts{1} = [-.1,1;.1,1];
% pfix = any fixed points to provide to the model (default = [])
% param(1) = base_spacing - edge length away from refined nodes (eg 0.1)
% param(2) = refine_ratio - relative refinement near points (eg. 10)
% param(3) = gradient     - transition slope of refinement (eg 0.1)
%
% Example: Three electrodes 2 on boundary and one internal
%  elec_pts = {[1,0],[0,1;sin(0.2),cos(0.2)],[0.5,0.5]};
%  fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.15,10,0.05] );
%
% Example:
%  n_elecs= 14; elec_width= 0.1; hw= elec_width/2;
%  th = linspace(0,2*pi,n_elecs+1); th(end)=[];
%  for i=1:n_elecs;
%     ti = th(i) + [hw;-hw];
%     elec_pts{i} = [sin(ti),cos(ti)];
%  end
%  fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.02] );
%
% See also: dm_2d_pt_elecs

% (C) 2009 Andy Adler. License: GPL version 2 or version 3
% $Id: dm_2d_circ_pt_elecs.m 6926 2024-05-31 15:34:13Z bgrychtol $


if ischar(elec_pts) && strcmp(elec_pts,'UNIT_TEST'); do_unit_test; return; end
if ischar(elec_pts) && strcmp(elec_pts,'CALIBRATE'); do_calibrate; return; end

copt.cache_obj = {elec_pts, spacing};
copt.fstr = 'dm_2d_circ_pt_elecs';
fmdl = eidors_cache(@do_circ_pt_elecs,{elec_pts, pfix, spacing}, copt);

function fmdl = do_circ_pt_elecs( elec_pts, pfix, spacing )
params.base_spacing = spacing(1);
params.refine_ratio = spacing(2);
params.gradient     = spacing(3);

bbox= [-1,-1;1,1];

eidors_cache('boost_priority',-4);
fmdl= dm_2d_pt_elecs( elec_pts, [], params, @circle, [-1,-1;1,1] );
eidors_cache('boost_priority',+4);

fmdl.name = sprintf('dm_2d_circ_pt_elec');

function d= circle(p,params);
  d = sqrt(sum(p.^2,2)) - 1; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% TESTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function do_unit_test
   elec_pts = {[1,0],[0,1;sin(0.2),cos(0.2)],[0.5,0.5]};
   fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.15,10,0.05] );
   subplot(221); show_fem(fmdl);
   for i=1:length(elec_pts)
      unit_test_cmp(sprintf('Elec%d(%d)',1,i), ...
         mean(fmdl.nodes(fmdl.electrode(i).nodes,:),1), mean(elec_pts{i},1), .01);
   end


% Example:
   n_elecs= 14; elec_width= 0.1; hw= elec_width/2;
   th = linspace(0,2*pi,n_elecs+1); th(end)=[];
   for i=1:n_elecs;
      ti = th(i) + [hw;-hw];
      elec_pts{i} = [sin(ti),cos(ti)];
   end
   fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [0.10,10,0.02] );
   subplot(222); show_fem(fmdl)
   for i=1:length(elec_pts)
      unit_test_cmp(sprintf('Elec%d(%d)',2,i), ...
         mean(fmdl.nodes(fmdl.electrode(i).nodes,:),1), mean(elec_pts{i},1), .01);
   end
   eidors_msg('CHECK FIGURE - should have 3 and 16 electrodes',0);

% find example models that work well with the values defined in mk_common_model
function do_calibrate

% meas = calc_range(0.05*([1,1.1,1.2,1.4,1.6,2.0,2.5,3,4,5,7,10,15,20]),0,1) 
  meas = calc_range(0.02,3,0.10),%- level 1
% meas = calc_range(0.03,10,0.05),%- level 3

% meas = calc_range(0.02,10,0.05),%- level 3
% meas = calc_range(0.02,20,0.05),%- level 3
 % As a start, set p2 = 1. This is uniform, and we can choose the
 % target refinement levels, and calculate the centre mesh density
 % OR setting p3=1 also gives uniform

function meas = calc_range(p1range,p2,p3)
  n_elecs= 8; elec_width= 0.1; hw= elec_width/2;
  th = linspace(0,2*pi,n_elecs+1); th(end)=[];
  for i=1:n_elecs;
     ti = th(i) + [hw;-hw];
     elec_pts{i} = [sin(ti),cos(ti)];
  end

  for i = 1:length(p1range)
     p1 = p1range(i);
     fmdl= dm_2d_circ_pt_elecs( elec_pts, [], [p1,p2,p3]);
     nd= fmdl.nodes; nn= size(nd,1);
     ne= size(fmdl.elems,1);
     nc= sum( nd(:,1).^2 + nd(:,2).^2 < 0.5^2);
     nl= sum( (nd(:,1)-1).^2 + nd(:,2).^2 < 0.2^2);
     meas(i,:) = [p1,p2,p3,[nn,ne,nc,nl]/1000];

  end;
show_fem(fmdl)

