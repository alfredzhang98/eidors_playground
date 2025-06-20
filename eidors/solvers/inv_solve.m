function img = inv_solve( inv_model, data1, data2)
% INV_SOLVE: calculate imag from an inv_model and data
% 
% inv_solve can be called as
%     img= inv_solve( inv_model, data1, data2)
%   if inv_model.reconst_type = 'difference'
% or
%     img= inv_solve( inv_model, data )
%   if inv_model.reconst_type = 'static'
%
%   if inv_model.reconst_to = 'nodes' then output
%      img.node_data has output data
%   else   reconst_to = 'elems' (DEFAULT) then output to
%      img.elem_data has output data
%
% in each case it will call the inv_model.solve
%
% data      is a measurement data structure
% inv_model is a inv_model structure
% img       is an image structure
%           or a vector of images is data1 or data2 are vectors
%
% For difference EIT:
% data1      => difference data at earlier time (ie homogeneous)
% data2      => difference data at later time   (ie inhomogeneous)
%
% data can be:
%   - an EIDORS data object
%
%   - an M x S matrix, where M is the total number
%         of measurements expected by inv_model
%
%   - an M x S matrix, where M is n_elec^2
%        when not using data from current injection
%        electrodes, it is common to be given a full
%        measurement set.  For example, 16 electrodes give
%        208 measures, but 256 measure sets are common.
%        Data will be selected based on fwd_model.meas_select.
%
% If S > 1 for both data1 and data2 then the matrix sizes must be equal
%
% Parameters:
%   inv_model.inv_solve.select_parameters: indices of parameters to return
%                         DEFAULT: return all paramteres
%  Scale solution (to correct for amplitude or other defects)
%   inv_model.inv_solve.scale_solution.offset
%   inv_model.inv_solve.scale_solution.scale
%  Disable solution error calculations
%   inv_model.inv_solve.calc_solution_error = 0
%
% img.time is frame sequence number. It is in seconds if
%   imdl.fwd_model.frame_rate is specified


% (C) 2005-2010 Andy Adler. License: GPL version 2 or version 3
% $Id: inv_solve.m 6665 2024-03-08 18:41:02Z aadler $

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

inv_model = prepare_model( inv_model );
opts = parse_parameters( inv_model );

print_hello(inv_model.solve);
try inv_model.solve = str2func(inv_model.solve); end

if opts.abs_solve
   if nargin~=2;
      error('only one data set is allowed for a static reconstruction');
   end
   
   fdata = filt_data(inv_model,data1);

   imgc= eidors_cache( inv_model.solve, {inv_model, fdata}, 'inv_solve');

else
   if nargin~=3;
      error('two data sets are required for a difference reconstruction');
   end

   % expand data sets if one is provided that is longer
   data_width= max(num_frames(data1), num_frames(data2));

   fdata1 = filt_data( inv_model, data1, data_width );
   fdata2 = filt_data( inv_model, data2, data_width );
   % TODO: Check if solver can handle being called with multiple data
   imgc= eidors_cache( inv_model.solve, {inv_model, fdata1, fdata2},'inv_solve');
   

end

check_parametrization_handling(inv_model,imgc);
     
img = eidors_obj('image', imgc );
% img = data_mapper(img,1);
if ~isfield(img,'current_params') % current_params is Bartek's units conversion control for img.elem_data.

  img.current_params = [];
end

% If we reconstruct with a different 'rec_model' then
%  put this into the img
if isfield(inv_model,'rec_model')
   img.fwd_model= inv_model.rec_model;
end


if opts.select_parameters;
   img.elem_data = img.elem_data( opts.select_parameters, :);
end

if isfield(img,'elem_data')
   img.time = 1:size(img.elem_data,2);
else
   img.time = 1:size(img.node_data,2);
end

if isfield(inv_model.fwd_model,'frame_rate');
   img.time = img.time / inv_model.fwd_model.frame_rate;
end

if ~opts.reconst_to_elems
  img.node_data= img.elem_data;
  img = rmfield(img,'elem_data');
end

% Scale if required
try; img.elem_data = opts.offset + opts.scale * img.elem_data; end
try; img.node_data = opts.offset + opts.scale * img.node_data; end

% MATLAB IS SUPER SUPER STUPID HERE. YOU CAN'T ASSIGN IF FIELDS ARE IN
%  A DIFFERENT ORDER. Example
%>> A.a= 1; A.b= 2; B.b= 3; B.a = 3; A(2) = B
%??? Subscripted assignment between dissimilar structures.
img.info.error = NaN;
img = orderfields(img);

% calculate residuals
try 
   do_calc = inv_model.inv_solve.calc_solution_error;
catch
   eidors_msg('inv_solve: not calculting solution residual',3);
   do_calc = false;
end
if ~do_calc;
   return;
end
eidors_msg('inv_solve: Calculating solution residual (inv_model.inv_solve.calc_solution_error = 0 to disable)',2);
try
   if opts.abs_solve
      img.info.error = calc_solution_error( imgc, inv_model, fdata);
   else
      img.info.error = calc_solution_error( imgc, inv_model, fdata1, fdata2 );
   end
   eidors_msg('inv_solve: Solution Error: %f', img.info.error,  2);
catch
   eidors_msg('inv_solve: Solution Error calculation failed.',2);
end

function print_hello(solver)
    if isa(solver,'function handle')
        solver = func2str(solver);
    end
    if strcmp(solver, 'eidors_default');
        solver = eidors_default('get','inv_solve');
    end
    eidors_msg('inv_solve: %s', solver,3);


function opts = parse_parameters( imdl )
   if  strcmp(imdl.reconst_type,'static') || ...
       strcmp(imdl.reconst_type,'absolute')
      opts.abs_solve = 1;
   elseif strcmp(imdl.reconst_type,'difference')
      opts.abs_solve = 0;
   else
      error('inv_model.reconst_type (%s) not understood', imdl.reconst_type); 
   end

   opts.select_parameters = [];
   try
      opts.select_parameters = imdl.inv_solve.select_parameters;
   end

   opts.reconst_to_elems = 1;
   try if strcmp( imdl.reconst_to, 'nodes' )
      opts.reconst_to_elems = 0;
   end; end
   
   opts.scale  = 1;
   try opts.scale = imdl.inv_solve.scale_solution.scale; end

   opts.offset = 0;
   try opts.offset = imdl.inv_solve.scale_solution.offset; end

function imdl = deprecate_imdl_parameters(imdl)
   if isfield(imdl, 'parameters')
      if isstruct(imdl.parameters)
%          disp(imdl)
%          disp(imdl.parameters)
%          warning('EIDORS:deprecatedParameters','INV_SOLVE inv_model.parameters.* are deprecated in favor of inv_model.inv_solve.* as of 30-Apr-2014.');

         if ~isfield(imdl, 'inv_solve')
            imdl.inv_solve = imdl.parameters;
         else % we merge
            % merge struct trick from:
            %  http://stackoverflow.com/questions/38645
            A = imdl.parameters;
            B = imdl.inv_solve;
            M = [fieldnames(A)' fieldnames(B)'; struct2cell(A)' struct2cell(B)'];
            try % assumes no collisions
               imdl.inv_solve=struct(M{:});
            catch % okay, collisions - do unique to resolve them
               [tmp, rows] = unique(M(1,:), 'last');
               M=M(:,rows);
               imdl.inv_solve=struct(M{:});
            end
         end
         imdl = rmfield(imdl, 'parameters');
         imdl.parameters = imdl.inv_solve; % backwards compatible!
      else
         error('unexpected inv_model.parameters where parameters is not a struct... i do not know what to do');
      end
   end

function mdl = prepare_model( mdl )
    fmdl = mdl.fwd_model;
    fmdl = mdl_normalize(fmdl,mdl_normalize(fmdl));
    if ~isfield(fmdl,'elems');
        return;
    end

    fmdl.elems  = double(fmdl.elems);
    fmdl.n_elem = size(fmdl.elems,1);
    fmdl.n_node = size(fmdl.nodes,1);
    if isfield(fmdl,'electrode');
        fmdl.n_elec = length(fmdl.electrode);
    else
        fmdl.n_elec = 0;
    end

    mdl.fwd_model= fmdl;
    if ~isfield(mdl,'reconst_type');
        mdl.reconst_type= 'difference';
    end

    % warn if we have deprecated inv_model.parameters laying about
    mdl = deprecate_imdl_parameters(mdl);

function check_parametrization_handling(inv_model,imgc)
if isfield(inv_model, 'jacobian_bkgnd') && ... 
    has_params(inv_model.jacobian_bkgnd) && ~has_params(imgc)
   if isa(inv_model.solve,'function_handle')
      solver = func2str(inv_model.solve);
   else
      solver = inv_model.solve;
   end
   if strcmp(solver,'eidors_default');
      solver = eidors_default('get','inv_solve');
   end
   warning('EIDORS:ParamsObliviousSolver',...
      ['The solver %s did not handle the parametrization data properly.\n'...
       'The results may be incorrect. Please check the code to verify.'], ...
       solver);
end
         
    
function b = has_params(s)
b = false;
if isstruct(s)
   b = any(ismember(fieldnames(s),supported_params));
end


% TODO: this code really needs to be cleaned, but not before eidors 3.4
function nf= num_frames(d0)
   if isnumeric( d0 )
      nf= size(d0,2);
   elseif d0(1).type == 'data';
      nf= size( horzcat( d0(:).meas ), 2);
   else
      error('Problem calculating number of frames. Expecting numeric or data object');
   end
   
% test for existance of meas_select and filter data
function d2= filt_data(inv_model, d0, data_width )
   if ~isnumeric( d0 )
       % we probably have a 'data' object

       d1 = [];
       for i=1:length(d0)
          if strcmp( d0(i).type, 'data' )
              d1 = [d1, d0(i).meas];
          else
              error('expecting an object of type data');
          end
       end

   else
      % we have a matrix of data. Hope for the best
      d1 = d0;
   end

   d1= double(d1); % ensure we can do math on our object

   if isfield(inv_model.fwd_model,'meas_select') && ...
     ~isempty(inv_model.fwd_model.meas_select)
      % we have a meas_select parameter that isn []

      meas_select= inv_model.fwd_model.meas_select;
      if     size(d1,1) == length(meas_select)
         d2= d1(meas_select,:);
      elseif size(d1,1) == sum(meas_select==1)
         d2= d1;
      else
         error('inconsistent difference data: (%d ~= %d). Maybe check fwd_model.meas_select',  ...
               size(d1,1), length(meas_select));
      end
   else
      d2= d1;
   end

   if nargin==3 % expand to data width
      d2_width= size(d2,2);
      if d2_width == data_width
         % ok
      elseif d2_width == 1
         d2= d2(:,ones(1,data_width));
      else
         error('inconsistent difference data: (%d ~= %d)',  ...
               d2_width, data_width);
      end
   end

% Test code
function do_unit_test
   k=0; N=5; nd = 5;

   imdl = mk_common_model('d2c2',16);
   imdl = select_imdl( imdl, {'Choose NF=2.5'});
   mvx = linspace(-0.8,0.8,nd);
   [vh,vi] = simulate_movement(mk_image(imdl), [mvx;0*mvx;0.05+0*mvx]);

   img= inv_solve(imdl,vh,vi); img.show_slices.img_cols = 5;
   k=k+1; subplot(N,1,k); show_slices(img);
   unit_test_cmp('inv_solve: 1a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 1b',  std(img.elem_data), 1.5e-2, 1e-2);

   vhm= eidors_obj('data','nn','meas',vh);
   img= inv_solve(imdl,vhm,vi);
   unit_test_cmp('inv_solve: 2a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 2b',  std(img.elem_data), 1.5e-2, 1e-2);

   img= inv_solve(imdl,vh*ones(1,nd),vi);
   unit_test_cmp('inv_solve: 3a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 3b',  std(img.elem_data), 1.5e-2, 1e-2);

   vim= eidors_obj('data','nn','meas',vi);
   img= inv_solve(imdl,vhm,vim);
   unit_test_cmp('inv_solve: 4a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 4b',  std(img.elem_data), 1.5e-2, 1e-2);

   vhm= eidors_obj('data','nn','meas',vh*ones(1,nd));
   img= inv_solve(imdl,vhm,vi);
   unit_test_cmp('inv_solve: 5a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 5b',  std(img.elem_data), 1.5e-2, 1e-2);

   vhm(1)= eidors_obj('data','nn','meas',vh*ones(1,2));
   vhm(2)= eidors_obj('data','nn','meas',vh*ones(1,nd-2));
   img= inv_solve(imdl,vhm,vi);
   unit_test_cmp('inv_solve: 6a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 6b',  std(img.elem_data), 1.5e-2, 1e-2);

   vim(1)= eidors_obj('data','nn','meas',vi(:,1:3));
   vim(2)= eidors_obj('data','nn','meas',vi(:,4:end));
   img= inv_solve(imdl,vhm,vim);
   unit_test_cmp('inv_solve: 7a', mean(img.elem_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 7b',  std(img.elem_data), 1.5e-2, 1e-2);

   im2 = imdl; im2.inv_solve.select_parameters = 1:5;
   img= inv_solve(im2,vh,vi);
   unit_test_cmp('inv_solve: 10', size(img.elem_data), [5,nd]);


   im2 = imdl;
   im2.inv_solve.scale_solution.offset = 1;
   im2.inv_solve.scale_solution.scale = 2;
   img= inv_solve(im2,vh,vi); img.show_slices.img_cols = 5;
   unit_test_cmp('inv_solve: 20a', mean(img.elem_data), 1.006, 2e-3);
   unit_test_cmp('inv_solve: 20b',  std(img.elem_data), 3e-2, 1e-2);
   im2.inv_solve.calc_solution_error = 1;
   img= inv_solve(im2,vh,vi);
   unit_test_cmp('inv_solve: 20e', mean(img.elem_data), 1.006, 2e-3);
   
   im2.inv_solve.scale_solution.offset = 0;
   d = interp_mesh( imdl.fwd_model); d= sqrt(sum(d.^2,2));
   im2.inv_solve.scale_solution.scale = spdiags(1-d,0,length(d),length(d));
   img= inv_solve(imdl,vh,vi); img.show_slices.img_cols = 5;
   k=k+1; subplot(N,1,k); show_slices(img);

   im2 = select_imdl(imdl, {'Nodal GN dif'} );
   img= inv_solve(im2,vh,vi);
   unit_test_cmp('inv_solve: 30a', mean(img.node_data), 3e-3, 1e-3);
   unit_test_cmp('inv_solve: 30b',  std(img.node_data), 1.5e-2, 1e-2);

   im2.inv_solve.scale_solution.offset = 1;
   im2.inv_solve.scale_solution.scale = 2;
   img= inv_solve(im2,vh,vi); img.show_slices.img_cols = 5;
   unit_test_cmp('inv_solve: 31a', mean(img.node_data), 1.006, 2e-3);
   unit_test_cmp('inv_solve: 31b',  std(img.node_data), 3e-2, 1.5e-2);

   im2 = select_imdl(imdl, {'Basic GN abs'} );
   im2.inv_solve.calc_solution_error = 0; % ensure accepted
   img= inv_solve(im2,vi(:,1));
   unit_test_cmp('inv_solve: 40a', mean(img.elem_data), 1.004, 1e-3);
   unit_test_cmp('inv_solve: 40b',  std(img.elem_data), 1.5e-2, 1e-2);
   im2.inv_solve.calc_solution_error = 1;
   img= inv_solve(im2,vi(:,1));
   unit_test_cmp('inv_solve: 40e', mean(img.elem_data), 1.004, 1e-3);


