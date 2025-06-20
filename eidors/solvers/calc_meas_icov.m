function meas_icov = calc_meas_icov( inv_model )
% meas_icov = calc_meas_icov( inv_model )
% CALC_MEAS_ICOV: calculate inverse covariance of measurements
%   The meas_icov is a matrix n_meas x n_meas of the
%     inverse of measurement covariances. Normally measurements
%     are assumed to be independant, so the meas_icov is 
%     a diagonal matrix of 1/var(meas)
%
% calc_meas_icov can be called as
%    meas_icov= calc_meas_icov( inv_model )
%
% meas_icov   is the calculated data prior
% inv_model    is an inv_model structure
%
% if:
%    inv_model.meas_icov    is a function
%          -> call it to calculate meas_icov
%    inv_model.meas_icov    is a matrix
%          -> return it as meas_icov
%    inv_model.meas_icov    does not exist
%          -> use I, or 1./homg (for normalized difference)

% (C) 2005 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_meas_icov.m 6718 2024-04-02 18:15:20Z aadler $

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

if isfield(inv_model,'meas_icov')
   if isnumeric(inv_model.meas_icov);
      meas_icov = inv_model.meas_icov;
   else
      try inv_model.meas_icov = str2func(inv_model.meas_icov); end
      meas_icov= eidors_cache( inv_model.meas_icov, {inv_model});
   end
else
   meas_icov= eidors_cache(@default_meas_icov,{inv_model} );
end

 
% Calculate a data prior for an assumption of uniform noise
% on each channel
% 
function meas_icov = default_meas_icov( inv_model )

   fwd_model= inv_model.fwd_model;
   fwd_model = mdl_normalize(fwd_model,mdl_normalize(fwd_model));

   n =  calc_n_meas( fwd_model );

   if ~mdl_normalize(fwd_model);
      meas_icov= speye( n );
   else
      homg_data=  solve_homg_image( fwd_model );
% if we normalize, then small data get increased
%   this means that noise on small data gets increased,
%    so the covariance is large when data are small
%   so the icov is small when data are small
% sig = k/h -> std = k/h -> 1/std = kh
      meas_icov = sparse(1:n, 1:n, ( homg_data.meas ).^2 );
   end

function n_meas = calc_n_meas( fwd_model )

   n_meas = 0;
   for i= 1:length(fwd_model.stimulation );
       n_meas = n_meas + size(fwd_model.stimulation(i).meas_pattern,1);
   end

% create homogeneous image + simulate data
function homg_data = solve_homg_image( fwd_mdl )
    n_elems= size( fwd_mdl.elems , 1);
    mat= ones( n_elems, 1);
    homg_img= eidors_obj('image', 'homogeneous image', ...
                         'elem_data', mat, 'fwd_model', fwd_mdl );
    homg_data=fwd_solve( homg_img);

function do_unit_test
   mdl= getfield( mk_common_model('a2c2',16), 'fwd_model');;
   run_dataprior_test( mdl );

   mdl = mdl_normalize(mdl,1); mdl.name = [mdl.name,':normalize'];
   run_dataprior_test( mdl );

   mdl= getfield( mk_common_model('n3r2',[16,2]), 'fwd_model');;
   run_dataprior_test( mdl );
   mdl = mdl_normalize(mdl,1); mdl.name = [mdl.name,':normalize'];
   run_dataprior_test( mdl );

% test dataprior
%     Difference dataprior should be 1
%     normalized difference dataprior should be homg_data.^2
function ok= run_dataprior_test( mdl )
    img= mk_image(mdl,1);
    homg_data=fwd_solve( img);

    DP= calc_meas_icov( img );

    % difference dataprior
    testvec= diag(DP);
    if mdl_normalize(mdl)
%  from calc_meas_icov, we have the following
%     meas_icov = sparse(1:n, 1:n, ( homg_data.meas ).^2 );
        testvec = homg_data.meas.^2 ./ diag(DP);
    end

    unit_test_cmp(['dataprior ',mdl.name],testvec,1,1e-10);

