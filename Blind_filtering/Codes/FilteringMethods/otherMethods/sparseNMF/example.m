% Example script for Sparse NMF with beta-divergence distortion function
% and L1 penalty on the activations.
%
% If you use this code, please cite:
% J. Le Roux, J. R. Hershey, F. Weninger, 
% "Sparse NMF ? half-baked or well done?," 
% MERL Technical Report, TR2015-023, March 2015
%   @TechRep{LeRoux2015mar,
%     author = {{Le Roux}, J. and Hershey, J. R. and Weninger, F.},
%     title = {Sparse {NMF} -?half-baked or well done?},
%     institution = {Mitsubishi Electric Research Labs (MERL)},
%     number = {TR2015-023},
%     address = {Cambridge, MA, USA},
%     month = mar,
%     year = 2015
%   }
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Copyright (C) 2015 Mitsubishi Electric Research Labs (Jonathan Le Roux,
%                                         Felix Weninger, John R. Hershey)
%   Apache 2.0  (http://www.apache.org/licenses/LICENSE-2.0) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% You need to provide a non-negative matrix v to be factorized.

params = struct;

% Objective function
params.cf = 'kl';   %  'is', 'kl', 'ed'; takes precedence over setting the beta value
  % alternately define: params.beta = 1;
params.sparsity = 5;

% Stopping criteria
params.max_iter = 100;
params.conv_eps = 1e-3;
% Display evolution of objective function
params.display   = 0;

% Random seed: any value over than 0 sets the seed to that value
params.random_seed = 1;

% Optional initial values for W 
%params.init_w
% Number of components: if init_w is set and r larger than the number of
% basis functions in init_w, the extra columns are randomly generated
params.r = 500;
% Optional initial values for H: if not set, randomly generated 
%params.init_h

% List of dimensions to update: if not set, update all dimensions.
%params.w_update_ind = true(r,1); % set to false(r,1) for supervised NMF
%params.h_update_ind = true(r,1);

[w, h, objective] = sparse_nmf(v, params);
