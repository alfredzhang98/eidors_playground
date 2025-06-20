function R_prior = calc_R_prior( inv_model, varargin )
% R = calc_R_prior( inv_model, varargin )
% CALC_R_PRIOR: calculate regularization matrix R
%   The image prior is matrix n_elem x ??? 
% 
% NOTE: This function is mostly just a hack. It's not obvious
%    how to to this correctly.
%
% calc_R_prior can be called as
%    R_prior= calc_R_prior( inv_model, ... )
%
% and will call the function inv_model.R_prior
% parameters to R_prior should be passed in the field
% inv_model.R_prior_function_name.parameters
% 
% If inv_model.R_prior is a matrix, calc_R_prior will return that matrix,
% possibly correcting for coarse2fine
%
% R_prior      calculated regularization prior R
% inv_model    is an inv_model structure

% (C) 2005-2008 Andy Adler. License: GPL version 2 or version 3
% $Id: calc_R_prior.m 7033 2024-11-29 00:24:25Z aadler $

if ischar(inv_model) && strcmp(inv_model,'UNIT_TEST'); do_unit_test; return; end

inv_model = rec_or_fwd_model( inv_model);


if isfield(inv_model,'R_prior')
   if isnumeric(inv_model.R_prior)
      R_prior = inv_model.R_prior;
   else
      R_prior= eidors_cache( inv_model.R_prior, inv_model );
   end
elseif isfield(inv_model,'RtR_prior')
   R_prior = eidors_cache(@calc_from_RtR_prior, inv_model,'calc_R_prior');
else
   error('calc_R_prior: neither R_prior or RtR_prior func provided');
end

if isfield(inv_model.fwd_model,'coarse2fine')
   c2f= inv_model.fwd_model.coarse2fine;
   if size(R_prior,1)==size(c2f,1)
%     we need to take into account coarse2fine - using a reasonable tol
      f2c = c2f';         
      R_prior = f2c*R_prior*c2f;
   end
end

function R_prior = calc_from_RtR_prior(inv_model)
   % The user has provided an RtR prior. We can use this to
   % get R =RtR^(1/2). Not that this is non unique
   if isnumeric(inv_model.RtR_prior)
      RtR_prior = inv_model.RtR_prior;
   else
      RtR_prior= eidors_cache( inv_model.RtR_prior, inv_model );
   end
   
   if exist('ldl') % not available for octave
      [L,D,P] = ldl(RtR_prior);
      R_prior = sqrt(D)*L'*P';
   else
      R_prior = sparse(size(RtR_prior,1),size(RtR_prior,2));
      [R,f,P] = chol(RtR_prior); % will often be sigular
      if f==0; f=size(RtR_prior,1); end % successful
      if f==1; f=size(R,1); end
      R_prior(1:f,:) = R*P';
   end

function inv_model = rec_or_fwd_model( inv_model);

   if isfield(inv_model,'rec_model');
      use_rec_model = 1;
      try if inv_model.prior_use_fwd_not_rec== 1;
         use_rec_model = 0;
      end; end

      if use_rec_model
          % copy the normalize flag from the fwd_model to prevent warnings
         inv_model.rec_model = mdl_normalize(inv_model.rec_model, ...
                                       mdl_normalize(inv_model.fwd_model));
         inv_model.fwd_model= inv_model.rec_model;
         inv_model= rmfield(inv_model,'rec_model');
      end
   end

function do_unit_test

   priors={'prior_tikhonov',
           'prior_laplace',
%          'prior_gaussian_HPF', % This is an R prior
           'prior_movement'};
   models= {'a2c0','f2c2','n3r2'};

   for j= 1:length(models);
      imdl = mk_common_model(models{j},16);
      for i = 1:length(priors);
         imdl.RtR_prior = feval(priors{i},imdl);

         R   = calc_R_prior(imdl);
         RtR = calc_RtR_prior(imdl);  
         if exist('OCTAVE_VERSION') && strcmp(priors{i},'prior_movement')
            eidors_msg('prior_movement test not valid in octave (need ldl)');
         else
            unit_test_cmp([models{j},':',priors{i}], R'*R, RtR, 1e-10);
         end
      end
   end


