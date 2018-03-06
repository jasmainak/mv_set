function [tr_lab,ts_lab,penval,vol,sz,dp]=ddt_core(xtrain,ytrain,xtest,max_cuts,mode);
%
% DDT_CORE: Compute DDT minimizing additive cost function using algorithm
% of Blanchard et al, 2004, 2005.
%
% USAGE: 
%
%	ddt_core(xtrain,ytrain,xtest,max_cuts,mode);
%	ddt_core(xtrain,ytrain,xtest,max_cuts,mode,'switch_name',switch_value,...);
%
% INPUT:
%
%   xtrain: d x n_train training patterns
%   ytrain: 1 x n_train training labels (0 through num_class - 1)
%   xtest: d x n_test testing patterns
%   max_cuts: 1 x d list of max number of splits for each dimension
%       if scalar, then replicated along each dimension
%   mode: number 1 through 3 characterizing learning problem
%       1: PE, CS, RL and NP/MM with one penalty for all classes
%       2: NP, MM with penalties for each class
%       3: MV, DL
%
% SWITCHES:
%
%   penalty: 
%     1: constant penalty of one per leaf
%     2: Rademacher upper bound
%     3: exact Rademacher
%     4: Rad. upper bound with volume term (mode = 3 only)
%
%   pen_wt: 
%     positive scalar that multiplies penalty [default=1]
%     for NP, MM, can be a vector. can also be vector of length
%     num_classes, specifying different penalty weights for each class
%
%   costs: 
%     1 x n_train vector of costs assigned to training examples
%     Used for CS, NP, MM, MV, DL, and RL problems
%
%   expermode:
%     0 -> verbose mode (display sizes of dictionary at different depths)
%     1 -> quiet mode
%     2 -> debug mode (depreciated)
%
% OUTPUT
%
%   tr_lab: 1 x n_train train labels 
%   ts_lab: 1 x n_test test labels 
%   err: empirical error or errors. this is a vector of length num_classes+1
%       when mode = 1 or 2, the last num_classes entries are the costs for
%           the respective classes, while the first entry is the cost of
%           the overall set estimate
%       when mode = 3, the first element is the volume and the second is
%           the cost of the set estimate
%   pen: same dimensions and meaning as error, except when mode=3 the first
%       entry is meaningless. note that the returned penalty values are NOT
%       scaled by the 'pen_wt' argument if one was given
%   sz: size
%   dp: depth
%
% COMPILING
%
%   Type "mex ddt_core.cpp" in the current directory to create the
%   executable if necessary.
%
%   If you have not run mex before, you will need to run "mex -setup" to
%   choose a compiler. Be sure to choose a compiler that can handle C++,
%   such as Microsoft Visual Studio.
%
% NOTES:
%
%   Wrapped for Matlab by Jason Laska
%
%   2005, Rice University
%