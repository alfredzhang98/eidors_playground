function out_img= show_slices( img, levels, vh )
% out_img = show_slices (img, levels ) show slices at levels of an
%             using a fast rendering algorithm
% img    = EIDORS image struct, or a array of structs
%          or output of CALC_SLICES (or equivalent multi-dimensional
%          picture array to be processed by mk_mosaic
% out_img= matrix in the current colormap which we can use image(out_img);
%
% PARAMETERS:
%  levels = any level defintion accepted by level_model_slice
%  vh     = [Lx2] matrix specifying the horizontal (column), vertical (row)
%           position of each slice in the output image
% 
% IMAGE PARAMETERS:
%   img.show_slices.levels (same as above);
%   img.show_slices.img_cols = number of columns in image
%   img.show_slices.sep = number of pixels between the individual images
%   img.show_slices.do_colourbar = draw a colourbar (via calc_colours)
%   img.calc_colours.npoints = pixel width/height to map to
%   img.get_img_data.frame_select = which frames of image to display
%   img.show_slices.axes_msm = use mdl_slice_mapper for x,y axes
%        (only works for single images and when level is specified on the 
%        image struct)
%        %% Example
%        img.fwd_model.mdl_slice_mapper = struct('level',[inf,0,inf], ...
%          'x_pts', linspace(-2,2,50), 'y_pts',linspace(2,10,100));
%        img.show_slices.axes_msm = true; show_slices(img);
%   img.show_slices.contour_levels = true => Do a contour on the image
%   img.show_slices.contour_levels = #    => Put this many contour lines
%   img.show_slices.contour_levels = vector => Put contours at these locations
%   img.show_slices.contour_properties => e.g. {'Color',[0,0,0],'LineWidth',2}
%   img.calc_slices.filter => e.g. conv2(ones(5),ones(5))/5^4   
%
% Show slices is roughly equivalent to:
%   rimg = calc_slices(img,levels); 
%   rimg = calc_colours(rimg,img); image(rimg);
%
% Deprecated interface:
%  if levels is [L x 5] then levels= [x,y,z,h,v] where,
%           x,y,z specify the axes intercepts, and 
%           h,v   specify the horizontal, vertical position
%                 of that slice in the output image
%
% See also: LEVEL_MODEL_SLICE, CALC_SLICES, MK_MOSAIC, CALC_COLOURS

% (C) 2005-2024 Andy Adler. License: GPL version 2 or version 3
% $Id: show_slices.m 7090 2024-12-20 20:25:43Z aadler $

if ischar(img) && strcmp(img,'UNIT_TEST'); do_unit_test; return; end

sep = 0;
try sep = img(1).show_slices.sep;
end


% TODO: because we wanted (back in 2005) to let show_slices
% handle lots of different scenarios (calling with images and
% without. It is now crufty ... and should be changed.
do_calc_slices = 0;
if isstruct(img)
   if strcmp(img(1).type,'image');
      do_calc_slices= 1;
   else
      error('show_slices only works with type image');
   end;
end

if nargin<=1;
   try   
       levels = img(1).show_slices.levels;
   catch
       levels = [];
   end
end

% Levels are set in calc_slices for 2D models. -- removed 2024-09-27 AA
%if isempty(levels) && do_calc_slices %&& size(img(1).fwd_model.nodes,2)==2
%   levels= [Inf,Inf,0];
%end

if nargin < 3
    vh = [];
end
if isempty(vh) && isnumeric(levels) && size(levels,2) == 5
   vh = levels(:,4:5);
   levels = levels(:,1:3);
   warning('EIDORS:DeprecatedInterface', ...
    'Specifing 5-element level is deprecated and will cause an error in a future release.');
end

if do_calc_slices
   rimg= calc_slices( img, levels );
else
   rimg= img;
   if isstruct(levels)
       img = levels;
       levels = [];
   end
end
% Eventualy we need to add a filter control to do this -aa'19
%% This is replaced with calc_slices.filter
%if 0
%   filt = ones(5);
%    filt = conv2(filt,filt);
%   rnimg = rimg;
%   rimg(isnan(rimg)) = 0;
%   rimg = conv2(rimg, filt,'same') / sum(filt(:));
%   rimg(isnan(rnimg)) = NaN;
%end

n_col = 0;
try  n_col = img(1).show_slices.img_cols;
end

r_img = mk_mosaic(rimg, sep, vh, n_col);

c_img = calc_colours( r_img, img);
out_img= reshape(c_img, size(r_img,1), size(r_img,2) ,[]);


axes_msm = false;
pts = [];
try axes_msm = img.show_slices.axes_msm; end
try pts = mdl_slice_mapper(img.fwd_model,'get_points'); end 

if axes_msm && size(rimg,3)==1 && ~isempty(pts)
   msm.x_pts = pts{1}; msm.y_pts = pts{2};
   image(msm.x_pts, msm.y_pts, out_img);
   set(gca,'Ydir','normal');
else
   image(out_img);
   axis image
end
axis off
axis equal
axis tight

do_colourbar = false;
try do_colourbar = img.show_slices.do_colourbar; end
if do_colourbar
    calc_colours( r_img, img, 1);
end

% Do a contour plot?
if isfield(img(1),'show_slices') && isfield(img(1).show_slices,'contour_levels');
   clevs = img.show_slices.contour_levels;
   if isfield(img.show_slices,'contour_properties');
      contour_properties = img.show_slices.contour_properties;
   else
      contour_properties = {'Color',[0.2,0.2,0.3]};
   end

   if ~axes_msm;
      msm.x_pts = 1:size(rimg,2);
      msm.y_pts = 1:size(rimg,1);
   end
   ish= ishold;
   if isnumeric(clevs)
      if ~ish; hold on; end 
      contour(msm.x_pts, msm.y_pts, rimg, clevs, contour_properties{:});
      if ~ish; hold off; end 
   elseif clevs % true but not numeric
      if ~ish; hold on; end 
      if exist('OCTAVE_VERSION');
      % Octave contour doesn't accept props 
         contour(msm.x_pts, msm.y_pts, rimg);
      else
         contour(msm.x_pts, msm.y_pts, rimg, contour_properties{:});
      end
      if ~ish; hold off; end 
   else
      error('img.show_slices.contour_levels parameter not understood');
   end
end

if nargout==0; clear('out_img'); end

function do_unit_test
   clf; sp=0;
   %1
   img=calc_jacobian_bkgnd(mk_common_model('a2t3',8)); 
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1);
   sp=sp+1;subplot(4,5,sp); show_slices(img) 
   %2
   img.calc_colours.npoints= 128;
   sp=sp+1;subplot(4,5,sp); show_slices(img) 
   %3
   img.calc_colours.npoints= 32;
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1:3);
   sp=sp+1;subplot(4,5,sp); show_slices(img) 

   %4
   img.show_slices.img_cols= 1;
   sp=sp+1;subplot(4,5,sp); show_slices(img) 

   %5
   imgn = rmfield(img,'elem_data');
   imgn.node_data=toeplitz(1:size(img.fwd_model.nodes,1),1);

   img.elem_data = img.elem_data(:,1);
   img.fwd_model.mdl_slice_mapper.npx = 10;
   img.fwd_model.mdl_slice_mapper.npy = 20;
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   sp=sp+1;subplot(4,5,sp); show_slices(img);

   %6
   img.elem_data = img.elem_data(:,1);
   img.fwd_model.mdl_slice_mapper = rmfield(img.fwd_model.mdl_slice_mapper, {'npx','npy'});
   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-100,100,20);
   img.fwd_model.mdl_slice_mapper.y_pts = linspace(-150,150,30);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   sp=sp+1;subplot(4,5,sp); show_slices(img);

   %7
   sp=sp+1;subplot(4,5,sp); show_slices(imgn) 

   %8
   imgn.fwd_model.mdl_slice_mapper.x_pts = linspace(-100,100,20);
   imgn.fwd_model.mdl_slice_mapper.y_pts = linspace(-150,150,30);
   imgn.fwd_model.mdl_slice_mapper.level = [inf,inf,0];
   sp=sp+1;subplot(4,5,sp); show_slices(imgn) 


% 3D images
   %9
   img=calc_jacobian_bkgnd(mk_common_model('n3r2',[16,2])); 
   img.calc_colours.npoints= 16;
   img.elem_data=toeplitz(1:size(img.fwd_model.elems,1),1);
   sp=sp+1;subplot(4,5,sp); show_slices(img,2) 

   %10
   img.elem_data=img.elem_data*[1:3];
   sp=sp+1;subplot(4,5,sp); show_slices(img,2) 

   %11
   img.elem_data=img.elem_data(:,1:2);
   sp=sp+1;subplot(4,5,sp); show_slices(img,[inf,inf,1;0,inf,inf;0,1,inf]);

   %12
   img.show_slices.sep = 5;
   img.fwd_model.mdl_slice_mapper.x_pts = linspace(-1,1,20);
   img.fwd_model.mdl_slice_mapper.y_pts = linspace(-1,1,30);
   img.fwd_model.mdl_slice_mapper.level = [inf,inf,0];

   sp=sp+1;subplot(4,5,sp); show_slices(img,2) 

   %13
   img.fwd_model = rmfield(img.fwd_model, 'mdl_slice_mapper');
   img.elem_data = img.elem_data(:,1);
   levels=[inf,inf,1,1,1;
           0,inf,inf,2,1;
           0,1,inf,3,1];
   sp=sp+1;subplot(4,5,sp); show_slices(img,levels) 

   %14
   levels=[inf,inf,1,1,1;
           0,inf,inf,2,2;
           0,1,inf,  1,3];
   sp=sp+1;subplot(4,5,sp); show_slices(img,levels) 

   %15
   img.elem_data = img.elem_data * [1,2];
   levels=[inf,inf,1,1,1;
           0,inf,inf,2,1;
           0,1,inf,  3,1];
   sp=sp+1;subplot(4,5,sp); show_slices(img,levels) 
   
   %16
   sp=sp+1;subplot(4,5,sp); show_slices(img,levels(:,1:3),levels(:,4:5)) 
   
   %17
   m = calc_slices(img,levels(:,1:3));
   sp=sp+1;subplot(4,5,sp); show_slices(m) 
   
   %18
   img.elem_data = img.elem_data(:,1);
   img.show_slices.contour_levels = true;
   sp=sp+1;subplot(4,5,sp); show_slices(img) 

   %19
   img.fwd_model.mdl_slice_mapper = struct('level',[inf,inf,1], ...
     'x_pts', linspace(-1,1,50), 'y_pts',linspace(-2,2,100));
   img.show_slices.axes_msm = true;
   img.show_slices.contour_properties = {'LineWidth',2};
   img.show_slices.contour_levels = 1:300;
   sp=sp+1;subplot(4,5,sp); show_slices(img) 

   %20
   l48 = linspace(0.5,48.5,24+1); 
   l32 =-linspace(0.5,32.5,16+1) + 20; 
   rmfn = @(xyz) vecnorm(xyz,2,2)>40;
   fmdl = mk_grid_model([],l48,l32,l48,rmfn);
   xyz = interp_mesh(fmdl); 
   img = mk_image(fmdl, sum(xyz,2));
   sp=sp+1;subplot(4,5,sp); show_slices(img,3) 
