function hh= eidors_colourbar(max_scale,ref_lev, cb_shrink_move, greyscale)
% EIDORS_COLOURBAR - create an eidors colourbar with scaling to image
% usage: eidors_colourbar( img )
%    show a colourbar on the current axis representing img
%
% usage: hh= eidors_colourbar( img )
%    return a handle to colourbar created 
%    handle can be modified: set(hh,'ylim', ...)
%
% usage: eidors_colourbar( img );
%   img.calc_colours.ref_level =  %  centre of the colour scale
%   img.calc_colours.clim      =  %  max diff from ref_level
%
% usage: eidors_colourbar( img ); % Set own tick positions
%   img.eidors_colourbar.tick_vals = [-20:20]/10;
%
% usage: eidors_colourbar( img ); % Set own tick positions
%   img.eidors_colourbar.tick_divisions = 5; % Extra tick divisions
%
% Optional parameter:
%    cb_shrink_move(1) = horizontal shrink (relative)
%    cb_shrink_move(2) = vertial shrink (relative)
%    cb_shrink_move(3) = horizontal move (absolute screen units)
%  img.calc_colours.cb_shrink_move = [1,2,0];
%
% Make sure you use 'axis tight' in order for cb_shrink_move to
%    give a correct looking image
%
% KNOWN ISSUE: if you use cb_shrink_move, then matlab will
%   forget the link between the figure and its colorbar. Future
%   plots in the same axis will continue to shrink. In general, the
%   axis will need to be cleared or reinitialized.
% EXAMPLE:
%   show_slices(img,2);
%    p = get(gca,'position') 
%   eidors_colourbar(img);
%    set(gca,'position',p); %%% Reset axes after colourbar and move
%
% The colorbars are removed with colorbar('delete')

% (C) 2005-2010 Andy Adler. License: GPL version 2 or version 3
% $Id: eidors_colourbar.m 6562 2023-02-16 13:03:39Z aadler $


if ischar(max_scale) && strcmp(max_scale,'UNIT_TEST'); do_unit_test; return; end

% if called as a simple eidors colourbar function
if isstruct(max_scale) && strcmp(max_scale.type,'image')
% OLD WAY WAS TO CALL calc_colours first
%   calc_colours(max_scale,[],1);
%   return
   img = max_scale;
   pp=get_colours( img );
   img_data= get_img_data( img );
   [scl_data, ref_lev, max_scale] = scale_for_display( img_data(:), pp);
   try if strcmp(pp.component,'imag');
       ref_lev = imag(ref_lev);
   end; end
   ref_lev = real(ref_lev);
end

   hh= colorbar;
   drawnow

   cbsm = [];
   try;  cbsm =  img.calc_colours.cb_shrink_move;
   end;

   if nargin >= 3; cbsm = cb_shrink_move;
   end

   axpos = get(gca,'Position');
   if ~isempty(cbsm)
      do_cb_shrink_move( hh, cbsm );
   end

   %FIXME = AA+CG 30/1/12
   if nargin <4; 
       greyscale=[]; 
   else
       warning(['eidors_colourbar: greyscale is an experimental feature'...
           'and will be re-implemented']);
   end

   % Stop scale from being too small
   if max_scale<abs(ref_lev)
      if max_scale < 1e-10; max_scale = 1e-10; end
   else
      if max_scale/abs(ref_lev) < 1e-4; max_scale = ref_lev*1e-4; end 
   end

   % Get colormap limits  and move bottom so we don't see the background colour 
   ylim = get(hh,'Ylim');
   cmsz = size(colormap,1);
   cbsz = ylim(2) - ylim(1);
   unit = cbsz/cmsz;
   ylim(1)= ylim(1)+unit;
   %FIXME = AA+CG 30/1/12
   if ~isempty(greyscale); ylim(1) = ylim(1) + unit ; end
   set(hh,'Ylim',ylim);

   c_ctr = mean(ylim);
   c_max = ylim(2) - c_ctr;


   try
      tick_vals = img.eidors_colourbar.tick_vals;
   catch
      try
         tick_div = img.eidors_colourbar.tick_divisions;
      catch
         tick_div = [];
      end
      tick_vals= get_tick_vals(max_scale, ref_lev, greyscale, tick_div);
   end

   % ref_lev goes to c_ctr. max_scale goes to c_max
%FIXME - need a switch to control use of max scale
   tick_locs = (tick_vals - ref_lev)/max_scale * c_max + c_ctr;
if isempty(greyscale) 
    tick_locs = (tick_vals - ref_lev)/max_scale * c_max + c_ctr;
else
    tick_locs = tick_vals*c_max*2/max_scale +2.5;
end
   % fix extremes so numbers are displayed
   tick_locs(abs(tick_locs-ylim(1))<eps) = ylim(1);
   tick_locs(abs(tick_locs-ylim(2))<eps) = ylim(2);

   set(hh,'YTick', tick_locs');
   set(hh,'YTickLabel', tick_vals');
% Code to right align YTickLabel
if 0
   lab = get(hh,'YTickLabel');
   for i=1:size(lab,1);
      nsp = sum(lab(i,:)==' ');
      lab(i,:) = lab(i,[end-(nsp-1:-1:0),1:(end-nsp)]);
   end
   set(hh,'YTickLabel',lab);
end

   if ~isempty(cbsm) % RESET OUR AXES
      if ~all(cbsm == [1,1,0]); 
         set(gca,'position',axpos);
      end
   end

   % Clear hh do it is not displayed, unless output requested
   if nargout ==0; clear hh; end

end

function do_cb_shrink_move( hh, cbsm )
   % make colourbar smaller and closer to axis

   axpos = get(gca,'Position');
   posn= get(hh,'Position');
   if ~all(cbsm == [1,1,0]); 
      posn = [posn(1) - cbsm(3), posn(2) + posn(4)*(1-cbsm(2))/2, ...
              posn(3) * cbsm(1), posn(4) * cbsm(2)];
      set(hh,'Position', posn );
      set(gca,'Position',axpos);
   end

end

% This function is copied from calc_colours. (AA: 23/3/2015)
% The plan is to eventually separate eidors colourbar from
%  calc_colours in a future release. This is a step on the
%  way to separating the code
function pp=get_colours( img );
   global eidors_colours;
   pp= eidors_colours;

   pp.component = 'real'; % will be overriden if in image

% override global if calc.colours specified
   try
% DAMN Matlab should have syntax for this loop
      fds= fieldnames(img(1).calc_colours);%assume all are the same!
      for fdn= fds(:)';
         fdn= fdn{1};
         pp.( fdn ) = img(1).calc_colours.(fdn);
      end
   end
end

function tick_vals= get_tick_vals(max_scale, ref_lev, greyscale, tick_div_in) 
% COMMENTS: AA - 4 apr 13
% If we want to ahve less aggressive rounding, we can do this
   F= 2;
   OrdOfMag = 10^floor(log10(max_scale*F))/F;

   scale1  = floor( max_scale / OrdOfMag + 2*eps );
% disp([scale1, OrdOfMag, max_scale]); => DEBUG
   if     (scale1/F >= 8);  fms = F*8;   tick_div=2; 
   elseif (scale1/F >= 6);  fms = F*6;   tick_div=2; 
   elseif (scale1/F >= 4);  fms = F*4;   tick_div=2;
   elseif (scale1/F >= 3);  fms = F*3;   tick_div=3;
   elseif (scale1/F >= 2);  fms = F*2;   tick_div=2;
   elseif (scale1/F >= 1.5);fms = F*1.5; tick_div=3;
   elseif (scale1/F >= 1);  fms = F*1;   tick_div=2;
   else   (scale1/F >= 0.5);fms = F*0.5; tick_div=2;
   end

   if ~isempty(tick_div_in); tick_div = tick_div_in; end

%  ticks = (1:tick_div)/tick_div;
%  ticks(end) = [];

   scale_r  = OrdOfMag * fms;

%  in order to make the labels clean, we round to a near level
   OrdOfMag = 10^floor(log10(max_scale));
   ref_r = OrdOfMag * round( ref_lev / OrdOfMag );
   
   %FIXME = AA+CG 30/1/12
if isempty(greyscale)
    %  Set to 2 x range, so out of scale ticks are thrown
    tick_vals = linspace(-2,2,tick_div*4+1)*scale_r + ref_r;
else
%     tick_vals = [0:0.2:1]*max_scale;
     tick_vals = [0:0.2:1]*scale_r;
end
end

% DEBUG CODE ATTEMPTING TO FIX CB
function debug_code
         a = get(hh);
         set(hh,'Position', posn );
         a = rmfield(a,'CurrentPoint');
         a = rmfield(a,'TightInset');
         a = rmfield(a,'BeingDeleted');
         a = rmfield(a,'Type');
         a.Position = posn;
         set(hh,a);
         op = get(hh,'OuterPosition') 
         set(hh,'Position', posn );
         op1= get(hh,'OuterPosition') 
         set(hh,'OuterPosition',op);
         op2= get(hh,'OuterPosition') 
         set(gca,'Position',axpos);
end

function do_unit_test
   imgno = 1;
   imdl = mk_common_model('n3r2',[16,2]);
   img = mk_image(imdl);
   img=rmfield(img,'elem_data');
   img.node_data(1:252)= (0:251)/100 - 1.05;
   if 0
   fprintf('CMAX = %4.2f CMIN = %4.2f\n', ...
       max(img.node_data), min(img.node_data) );
   end


   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2); eidors_colourbar(img);

   % ref_level is not displayed, but will have its effect
   img.calc_colours.ref_level = -0.0234;
   subplot(3,3,imgno); imgno=imgno+1;
   img2= img;
   img2.eidors_colourbar.tick_divisions = 4;
   show_slices(img2,2); eidors_colourbar(img2);

   img.calc_colours.cmax = 1;
   subplot(3,3,imgno); imgno=imgno+1;
   img2= img;
   img2.eidors_colourbar.tick_vals = [-10:10]/5;
   show_slices(img2,2); eidors_colourbar(img2);

   MV =-0.05;
   img.calc_colours.cb_shrink_move = [.5,.8,MV];
   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2);
   eidors_colourbar(img);

   subplot(3,3,imgno); imgno=imgno+1;
   show_fem(img,1);

   subplot(3,3,imgno); imgno=imgno+1;
   show_slices(img,2);
    p = get(gca,'position');
   eidors_colourbar(img);
    set(gca,'position',p);

   show_slices(img,2);
    eidors_colourbar(img);

   subplot(3,3,imgno); imgno=imgno+1;
   img.node_data = 26e1*abs(img.node_data);
   show_slices(img,2);
   eidors_colourbar(img);

   subplot(3,3,imgno); imgno=imgno+1;
   img.node_data(1:252)= abs( (0:251)/100 - 1.05 );
   show_slices(img,2);
   eidors_colourbar(img);


   subplot(3,3,imgno); imgno=imgno+1;
   img.calc_colours.ref_level = 0.5;
   img.calc_colours.clim      = 0.5;
   show_slices(img,2);
   eidors_colourbar(img);

end
