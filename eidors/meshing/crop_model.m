function [fmdl,c2f_idx]= crop_model( axis_handle, fcn_handle );
% CROP_MODEL: Crop away parts of a fem model
%
% USAGE #1: crop display to show model internals
%   crop_model( axis_handle, fcn_handle );
%
%   fcn_handle ==1 where model is *kept*
% 
%   if axis_handle==[], then use the current axis
%   examples:
%     crop_model([],  inline('z==3','x','y','z'))
%     crop_model(gca, inline('x+y>0','x','y','z'))
%     crop_model([],  @(x,y,z) z<0);
%   if the fcn_handle is a string, a function with x,y,z is created
%     crop_model(gca, 'x+y>0') % same as previous
%
% USAGE #2: crop fwd_model to create new fwd_model
%   fmdl_new= crop_model( fwd_model, fcn_handle );
% 
%   example:
%   fmdl2= crop_model(fmdl1, @(x,y,z) x+y>0);
%
%  with two parameters output
% [fmdl,c2f_idx]= crop_model( axis_handle, fcn_handle );
%     c2f_idx maps each elemen in fmdl_new to fwd_model
%
% USAGE #3: crop img to create new img (preserve elem_data)
%   img2= crop_model(img1, @(x,y,z) x+y>0);

% (C) 2006-2008 Andy Adler. License: GPL version 2 or version 3
% $Id: crop_model.m 6811 2024-04-24 20:20:11Z bgrychtol $

if ischar(axis_handle) && strcmp(axis_handle,'UNIT_TEST'); do_unit_test; return; end

% TODO (update 2 apr 2012):
%  - make crop_model work for 2D fems

if isstr(fcn_handle)
  fcn_handle = inline(fcn_handle,'x','y','z');
end

type= isfield(axis_handle, 'type');
if type; type = axis_handle.type; end

if isempty(axis_handle)
   axis_handle= gca;
   crop_graphics_model(axis_handle, fcn_handle);
elseif ishandle( axis_handle )
   crop_graphics_model(axis_handle, fcn_handle);
elseif strcmp(type, 'fwd_model');
   [fmdl,c2f_idx]= crop_fwd_model(axis_handle, fcn_handle);
elseif strcmp(type, 'image');
   [fmdl_,c2f_idx]= crop_fwd_model(axis_handle.fwd_model, fcn_handle);
   fmdl = axis_handle; % input parameter
   fmdl.fwd_model = fmdl_;
   fmdl.elem_data = fmdl.elem_data(c2f_idx,:);
%  keyboard
else
   error('Not sure how to process first parameter');
end

% CROP GRAPHICS
function crop_graphics_model(axis_handle, fcn_handle);
   kk= get(axis_handle,'Children');
   % only the FEM frame
   %k=kk( find( kk== min(kk) ));

   for k= kk(:)'
      ktype = get(k,'Type');
      switch ktype; case {'patch','line'}
          crop_patch(fcn_handle, k)
      end % ignore non-patch
   end

function crop_patch(fcn_handle, k)
   x= get(k,'XData');
   y= get(k,'YData');
   try
      z= get(k,'ZData');
   catch
      z= [];
   end
   idx= ~all( feval(fcn_handle,x,y,z) );

   set_= {'Xdata', x(:,idx), 'Ydata', y(:,idx)};
   if ~isempty(z)
      set_{end+1}= 'Zdata';
      % Can't assign to empty array, so indirect
      dd = z(:,idx); 
      set_{end+1}= dd;
   end

   try
      c= get(k,'CData');
% Matlab doesn't document the shape of Cdata!
      if size(c,2)==1; c=c'; end
      if ~isempty(c);
          c = c(:,idx);
      end
%  This should work, but weird errors
%     set_= horzcat(set_,{'Cdata',c(:,idx)});
      set_{end+1}= 'Cdata';
      % Can't assign to empty array, so indirect
      set_{end+1}= c;
   end
   set(k, set_{:});

% CROP fwd_model
function [fmdl1,c2f_idx]= crop_fwd_model(fmdl0, fcn_handle);
   fmdl1= fmdl0;

% Find nodes to remove
   nodes= fmdl0.nodes;
   [n,d]= size(nodes);
   n2xyz= eye(d,3); 
   xyz= nodes*n2xyz;
% Keep these nodes
   idx0=  all( feval(fcn_handle,xyz(:,1), ...
                                xyz(:,2), ...
                                xyz(:,3)),2);
% remove these nodes
   [fmdl1, c2f_idx] = remove_nodes(fmdl0, idx0, 'keep');


function do_unit_test

   imdl = mk_common_model('a2c0',8); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<0','x','y','z'));
   unit_test_cmp('crop_model-a2c0-01',length(fmdl.electrode),8);
   unit_test_cmp('crop_model-a2c0-02',size(fmdl.elems),[32,3]);
   unit_test_cmp('crop_model-a2c0-03',size(fmdl.nodes),[25,2]);

   fmdl = crop_model(fmdl,'x<0'); % verify it's same
   unit_test_cmp('crop_model-str-a2c0-01',length(fmdl.electrode),8);
   unit_test_cmp('crop_model-str-a2c0-02',size(fmdl.elems),[32,3]);
   unit_test_cmp('crop_model-str-a2c0-03',size(fmdl.nodes),[25,2]);

   imdl = mk_common_model('n3r2',[16,2]); fmdl= imdl.fwd_model;
   fmdl = crop_model(fmdl,inline('x<0','x','y','z'));
   unit_test_cmp('crop_model-n3r2-01',length(fmdl.electrode),32);
   unit_test_cmp('crop_model-n3r2-02',size(fmdl.elems),[342,4]);
   unit_test_cmp('crop_model-n3r2-03',size(fmdl.nodes),[128,3]);
   fmdl = crop_model(fmdl,inline('z<2','x','y','z'));
   unit_test_cmp('crop_model-n3r2-04',length(fmdl.electrode),32);
   unit_test_cmp('crop_model-n3r2-05',size(fmdl.elems),[114,4]);
   unit_test_cmp('crop_model-n3r2-06',size(fmdl.nodes),[64,3]);




   subplot(331)
   imdl = mk_common_model('n3r2',[16,2]); fmdl= imdl.fwd_model;
   show_fem(fmdl);
   subplot(332)
   show_fem(fmdl); hh= gca; 
   subplot(333)
   show_fem(fmdl);
   crop_model([],inline('z<2','x','y','z'));
   crop_model(hh,inline('x<0','x','y','z'));

   subplot(334)
   fmdlc = fmdl;
   fmdlc = crop_model(fmdlc,inline('x<0','x','y','z'));
   show_fem(fmdlc);

   subplot(335)
   img = mk_image(fmdl,1); 
   load('datareal.mat','A'); img.elem_data(A)= 1.1;
   imgs =  crop_model(img,'y-z<-2.5'); % warning about # 11
   show_fem( imgs );
   unit_test_cmp('crop img',find( imgs.elem_data>1),[476;479; 482; 485])

   subplot(336)
   img = mk_image(fmdl,1); 
   load('datareal.mat','A'); img.elem_data(A)= 1.1;
   imgs =  crop_model(img,@(x,y,z) y-z<-2.5);
   show_fem( imgs );
   unit_test_cmp('crop img',find( imgs.elem_data>1),[476;479; 482; 485])



   subplot(337)
   imdl = mk_common_model('d2c2');
   fmdl= imdl.fwd_model;
   img = mk_image(imdl); img.elem_data(1:16) = 1.1;
   show_fem(fmdl);
   subplot(338)
   show_fem(img);
   try
   crop_model([],@(x,y,z) y<0)
   catch
   title('expected fail');
   end

   subplot(339)
   show_fem(fmdl); hh= gca; 
   crop_model(hh,inline('x<0','x','y','z'));


