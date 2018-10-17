function [labels,stats] = ddt_learn(xtrain,ytrain,xtest,params)
%
% DDT_LEARN     Set estimation with dyadic decision trees by voting over
%       multiple shifts and random forets
%
% USAGE
%
%   [labels,stats] = ddt_learn(xtrain,ytrain,xtest,params)
%
% INPUT
%
%   xtrain: d x n_train training data
%   ytrain: For supervised problems, 1 x n_train training labels, 0 through num_classes - 1
%           For unsupervised problems (MV and DL), input 1 x n_train vector of zeros
%   xtest: d x n_test test data
%   params: struct specifying parameters of learning algorithm
%
% Fields of params are as follows:
%
%   Fields applicable to all methods [default in brackets]
%   The first two must be specified
%
%     max_cuts: 
%       scalar or 1 x d vector; max number of cuts per dyadic rectangle per 
%         dimension
%       hint: if the dictionary size at the maximum depth (see display_dict 
%         parameter below) is close to the training sample size, your 
%         max_cuts is probably sufficient.
%     prob: string specifying learning problem. possible values are
%       'PE': minimum probability of error classification (supervised)
%       'CS': cost sensitive classification (supervised)
%       'NP': Neyman-Pearson classification (supervised)
%       'MM': min-max classification (supervised)
%       'MV': minimum volume sets (unsupervised)
%       'DL': density level sets (unsupervised)
%       'RL': regression level sets (supervised)
%     width_frac: [.99]
%       fraction of unit cube occupied by data, per feature
%     scale_meth: ['linear']
%       method for shrinking data down to unit cube
%       'linear': scale coordinates linearly
%       'quantile': space quantiles uniformly along each coordinate
%     n_forest: [1]
%       number of trees in random forest -- should only be > 1 when the
%       running time is too large otherwise.
%     rf_dim: [d]
%       number of features per tree in random forest -- should be as large
%       as possible within computational ability
%     n_shift: [1]
%       number of shifts to vote over, per tree of random forest -- in
%       general, the more shifts, the better the estimate
%     shift_frac: [0]
%       fraction of unit square to shift over in each coordinate;
%       set to 0 for unshifted DDTs
%     penalty: ['rad'] (note: 'rad' and 'con' seem to be the best performers)
%       'rad' = exact Rademacher average and assuming a uniform prior over
%           partitions
%       'radub' = upper bound on Rademacher average which is tight up to a
%           factor of 1/sqrt(2); useful if 'rad' is overfitting
%       'radvol' = Rademacher with volume term (MV/DL only). it was
%           experimentally observed (when d=2) that adding the volume to the 
%           empirical mass in the 'radub' penalty eliminated there being tiny 
%           islands of mass around isolated points. however, i've also
%           observed it causing this exact effect, so use with caution
%       'con' = constant penalty of 1/n_train per cell
%     pen_wt: [1] 
%       multiplicative factor applied to penalty; can be a vector of weights 
%       for NP/MM unless using constant penalty
%     display_dict: [0] flag for displaying dictionary sizes
%       1: diplay dictionary sizes for each depth; note: to get the
%           results to display "real time" in matlab, you need to call your
%           top level function as a "shortcut." you can create a shortcut
%           by dragging the script name from the command history to the
%           shortcut toolbar. also note that for NP/MM/MV, which require
%           iterative searches for the right lagrange multiplier, only the
%           dictionary sizes for the first iteration are displayed, since
%           they're all the same
%       0: quiet (display nothing)
%     display_est: [0] boolean flag
%       1: if d=2, graphically display set estimate
%       0: don't display set estimate. 
%     display_vol: [0 unless n_shift == 1] 
%       1: compute and display volume 
%       0: don't
%          note: if n_shift==1, the volume is exact as returned by
%          ddt_core. otherwise, use uniform grid to approximate.
%     n_grid: [256,40,16,10,6,5,4,3 for d=2:9 and 2 for d>9]
%       grid size (in each dimension) for plotting set estimate in 2d 
%          note: the larger the grid size the slower the code will run 
%     display_iter: [0] boolean flag (NP/MM/MV only)
%       1: display iterations of lagrange mult. search
%       0: don't
%
%   Additional fields for specific learning problems
%
%   Additional fields for CS
%
%     costs: (required)
%       1 x n_train vector of costs associated with each sample
%
%   Additional fields for iterative methods (NP, MM, MV)
%     max_iter: [10]
%       max number of iterations in Lagrange mult search
%     shift_inside: [1]
%         1: loop over shifts inside the loop over the lagrange multiplier
%         (recommended -- provides increased granularity in multiplier
%         search)
%         0: loop over shifts outside the loop over the lagrange multiplier
%
%   Additional fields for NP
%
%     alpha: [0.1] 
%       false alarm constraint
%     nu: [0]
%       damping factor on class 0 penalty, not to be confused with pen_wt
%     reweight_meth [2]
%       method for estimating/penalizing error:
%         1 -> treat as cost-sensitive problem and tune the cost parameter
%         2 -> bound class conditional errors separately and minimize
%           weighted sum of bounds
%
%   Additional fields for MM
%     rho: [1]
%       aim to have false positive rate = rho * false negative rate
%     reweight_meth [2]
%       method for estimating/penalizing error:
%         1 -> treat as cost-sensitive problem and tune the cost parameter
%         2 -> bound class conditional errors separately and minimize
%           weighted sum of bounds
%
%   Additional fields for MV
%
%     alpha: [0.9]
%       mass constraint
%     nu: [0]
%       damping factor, not to be confused with pen_wt
%
%   Additional fields for DL
%
%     level: (required)
%       targeted density level
%
%   Additional fields for RL
%
%     level: (required)
%       targeted regression function level
%
% OUTPUT:
% 
%   labels = 1 x n_test vector of labels
%   stats = struct whose fields are various statistics of the algorithm.
%       (currently this output has only been partially implemented)
%       The fields are as follows:
%           e_time: elapsed time 
%           emperr: empirical (training) error(s) (vector for NP, MM)
%           bound:  value of error bound (PE, CS, RL only)
%           vol:    volume of set estimate (MV,DL only)
%           sh:     shifts used
%           feat:   feature sets used for random forests
%           depth:  depths of trees; n_shift x n_forest x max_iter
%           size:   sizes of trees; n_shift x n_forest x max_iter
%           
%
% NOTES:
%
%   - multiple classes only supported for PE, CS, and MM
%   - the hard work is done by ddt_core.cpp, which implements the algorithm
%       found in
%           G. Blanchard, C. Schäfer, Y. Rozenholc, ``Oracle bounds and exact
%           algorithm for dyadic classification trees,'' In Proceedings of the 
%           17th. Conference on Learning Theory (COLT 2004). Springer 
%           Lecture Notes in Artificial Intelligence (3120), 378-392, 2004. 
%       type "help ddt_core" for more info
%       mex interface for ddt_core implemented by Jason Laska, Rice University
%   - currently the outputting of the learned tree(s) is not supported
%
% AUTHOR
%
% Clayton Scott, August 2005, Rice University

t0 = clock;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize some variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[d,n_train]=size(xtrain);
[dtest,n_test]=size(xtest);
num_classes = max(ytrain)+1;
if num_classes > 2
    error('voting not yet implemented for mult. classes -- see end of ddt_learn.m')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set defaults, check for errors, preprocess, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if d ~=dtest
    error('Training and test data have different dimensions')
end
if unique(ytrain) ~= [0:num_classes-1]
    error('training labels must be in the range 0 ... num_classes-1')
end
if ~isfield(params,'max_cuts')
    error('Specify params.max_cuts')
else
    max_cuts = params.max_cuts;
    if length(max_cuts) ~= d
        if length(max_cuts) == 1
            max_cuts=ones(1,d)*max_cuts;
        else
            error('params.max_cuts must be row vector of length 1 or number of dims')
        end
    end
    if (sum(max_cuts ~= round(max_cuts))>0) | (sum(max_cuts < 1)>0)
        error('elements of params.max_cuts must be positive integers')
    end
end
if ~isfield(params,'prob')
    error('Specify set estimation problem by setting params.prop -- type "help ddt" for help')
end
switch params.prob
    case 'PE'
        prob=1;
    case 'CS'
        prob=2;
    case 'NP'
        prob=3;
    case 'MM'
        prob=4;
    case 'MV'
        prob=5;
    case 'DL'
        prob=6;
    case 'RL'
        prob=7;
end
if num_classes ~= 2
    if (prob == 4 | prob == 7)
        error('Only two classes allowed for this learning problem');
    end
end
if ~isfield(params,'penalty') | strcmp(params.penalty,'rad')
    pen_num=3;  % default: Rademacher upper bound
elseif strcmp(params.penalty,'radub')
    pen_num=2;  % exact Rademacher
elseif strcmp(params.penalty,'radvol')
    pen_num=4;  % rad. vol
elseif strcmp(params.penalty,'con')
    pen_num=1;  % constant
else
    error('penalty not recognized')
end
if ~isfield(params,'pen_wt')
    if (prob ==3 | prob == 4)
        pen_wt = ones(1,num_classes);
    else
        pen_wt=1;
    end
    if pen_num==1
        pen_wt=pen_wt/n_train; % default for constant
    end
else
    if (prob ==3 | prob == 4) & length(params.pen_wt)==1
        pen_wt=ones(1,num_classes)*params.pen_wt;
    else
    	pen_wt=params.pen_wt;
    end
end
if strcmp(params.prob,'NP') | strcmp(params.prob,'MM')
    if ~isfield(params,'reweight_meth')
        if pen_num ~= 1
            reweight_meth = 2;  % bounds class cond. errors separately
        else
            reweight_meth = 1;
        end
    else
        reweight_meth = params.reweight_meth;
    end
    if reweight_meth == 2 & pen_num ==1
        error('Cannot use reweight_meth = 1 when pen_num = 1')
    end
end
% set penalty weight
if ~isfield(params,'pen_wt')
    if (prob ==3 | prob == 4) & reweight_meth == 2
        pen_wt = ones(1,num_classes);
    else
        pen_wt=1;
    end
else
    if (prob ==3 | prob == 4) & reweight_meth == 2 & length(params.pen_wt)==1
        pen_wt=ones(1,num_classes)*params.pen_wt;
    else
    	pen_wt=params.pen_wt;
    end
end
if pen_num==1
    pen_wt=pen_wt/n_train; % constant pen per cell is 1/n_train
end
if (prob ==3 | prob == 4) & reweight_meth == 2 
    if (length(pen_wt) ~= num_classes)
        error('pen_wt must be a vector of length num_classes')
    end
else    % for all other probs, pen_wt should be scalar
    if length(pen_wt) ~= 1
        error('pen_wt must be a scalar')
    end
end
% initialize some summary stats
switch params.prob
    case {'NP','MM'}
        emperr=zeros(1,num_classes);
        bound=zeros(1,num_classes);
    case {'MV','DL'}
        emperr=zeros(1,2);
        bound=0;            % not computed
    otherwise
        emperr=0;
        bound=0;
end
if ~isfield(params,'width_frac')
    wf=0.99;     % default
else
    wf=params.width_frac;
end
if ~isfield(params,'scale_meth')
    scale_meth='linear';    % default
else
    scale_meth=params.scale_meth;
end
if ~isfield(params,'delta')     % bound confidence
    delta=0.05;     
else
    delta=params.delta;
end

% ensemble params
if ~isfield(params,'rf_dim') | params.rf_dim == 0
    rf_dim=d;    
else
    rf_dim = params.rf_dim;
end
if ~isfield(params,'n_forest')
    n_forest = 1;   % default: no random forest
else
    n_forest = params.n_forest;
end
if n_forest*rf_dim < d
    warning('Not enough trees in random forest to cover all dimensions')
end
if n_forest > 1 & rf_dim == d
    error('Set params.rf_dim < d if params.n_forest > 1')
end
if ~isfield(params,'n_shift')
    n_shift = 1;     % default: don't vote over shifts
else
    n_shift = params.n_shift;
end
if ~isfield(params,'shift_frac')
    shift_frac = 0.5;
else
    shift_frac = params.shift_frac;
end
if n_shift ==1
    shift_frac = 0;  % if only considering one shift, make it unshifted
end
if shift_frac > 1 | shift_frac < 0
    error('params.shift_frac should be between 0 and 1');
end
if shift_frac == 0
    n_shift = 1; % larger n_shift is a waste of time
end
if ~isfield(params,'shift_type')
    shift_type = 'rand';
else
    shift_type = params.shift_type;
end

% params for specific probs
if strcmp(params.prob,'PE')
    costs = ones(1,n_train)/n_train;
end
if strcmp(params.prob,'CS') 
    if ~isfield(params,'costs')
        error('Set params.costs for cost sensitive classification')   
    else
        costs=params.costs;
        if length(find(costs < 0)) > 0
            error('costs must be non-negative')
        end
        costs=costs/sum(costs);
    end
end
if strcmp(params.prob,'NP') 
    if ~isfield(params,'alpha')
        warning('params.alpha set to 0.1 Neyman-Pearson classification')
        alpha=0.1;
    else
        alpha=params.alpha;
    end
end
if strcmp(params.prob,'NP') | strcmp(params.prob,'MV')
    if ~isfield(params,'nu')
        nu = 0;
    else
        nu = params.nu;
        if (nu < -1) | (nu > 1)
            error('nu should be between -1 and 1')
        end        
    end
end
if strcmp(params.prob,'NP') | strcmp(params.prob,'MV') | strcmp(params.prob,'MM')
    if ~isfield(params,'max_iter')
        max_iter = 10;
    else
        max_iter = params.max_iter;
    end
    if ~isfield(params,'tol')
        tol=0.0001;
    else
        tol = params.tol;
    end
    if ~isfield(params,'shift_inside')
        shift_inside = 1;
    else
        shift_inside = params.shift_inside;
    end
    if shift_inside == 1
        max_outer = max_iter;
        max_inner = 1;
    else
        max_outer = 1;
        max_inner = max_iter;
    end
else
    shift_inside = 0;
    max_outer = 1;
    max_inner = 1;
end
if strcmp(params.prob,'MM') 
    if ~isfield(params,'rho')
        rho=1;
    else
        rho=params.rho;
    end
    if rho <=0
        error('params.rho must be > 0')
    end
end
if strcmp(params.prob,'MV') 
    if ~isfield(params,'alpha')
        warning('params.alpha set to 0.9 for minimum volume set estimation')
        alpha=0.9;
    else
        alpha=params.alpha;
    end
end
if strcmp(params.prob,'DL') & ~isfield(params,'level')
    error('Set params.level density level set estimation')
end
if strcmp(params.prob,'RL') 
    if ~isfield(params,'level')
        error('Set params.level for regression level set estimation')
    end
    level=params.level;
    costs = abs(ytrain-level);  % will solve via cost sensitive classification
    costs = costs/sum(costs);
end
if strcmp(params.prob,'MV') | strcmp(params.prob,'DL')
    ytrain=zeros(1,n_train);
end
if strcmp(params.prob,'NP') | strcmp(params.prob,'MM')
    i0 = find(~ytrain); % useful variable
    n0 = length(i0);
    i1 = find(ytrain);
    n1 = length(i1);
end
if strcmp(params.prob,'NP')
    if nu==0
        params.tol=1/n0;
    end
end
if strcmp(params.prob,'MM')
    for i=1:num_classes
        ind{i}=find(ytrain==i-1);
        nsamp{i}=length(ind{i});
    end
end
if strcmp(params.prob,'MV')
    if nu == 0
        params.tol = 1/n_train;
    end
end

% display variables
if ~isfield(params,'display_est') | (num_classes > 2) | (d > 3)
    display_est=0;
else
    display_est=params.display_est;
end
if (~strcmp(params.prob,'MV') & ~strcmp(params.prob,'DL')) | ~isfield(params,'display_vol')
    if n_shift > 1
        display_vol=0;
    else
        display_vol=1;
    end
else
    display_vol = params.display_vol;
end
if display_est == 1 | display_vol ==1
    if ~isfield(params,'n_grid')
        switch d
            case 2
                n_grid=256;
            case 3
                n_grid=40;
            case 4
                n_grid=16;
            case 5
                n_grid=10;
            case 6
                n_grid=6;
            case 7
                n_grid=5;
            case 8
                n_grid=4;
            case 9
                n_grid=3;
            otherwise
                n_grid=2;
        end
    else
        n_grid=params.n_grid;
    end
else
    n_grid=0;
end
if ~isfield(params,'display_iter')
    display_iter=0;
else
    display_iter=params.display_iter;
end
if ~isfield(params,'display_dict')
    display_dict=0;    % quiet
else
    display_dict=params.display_dict;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normalize data to unit square
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch scale_meth
    case 'linear'
        %x=[xtrain];
        x=[xtrain xtest];     % "semi-supervised" normalization     
        max_vec=max(x,[],2);
        min_vec=min(x,[],2);
        width=max_vec-min_vec;
        xtrain=wf*(xtrain-repmat(min_vec,[1,n_train]))./repmat(width,[1,n_train]) + (1-wf)/2;
        xtest=wf*(xtest-repmat(min_vec,[1,n_test]))./repmat(width,[1,n_test]) + (1-wf)/2;
    case 'quantile'
        x=[xtrain xtest];     % "semi-supervised" normalization     
        tmp=linspace((1-wf)/2,1-(1-wf)/2,2*(n_train+n_test)+1);
        new=tmp(2:2:end);
        for i=1:d
            [xsort,ind]=sort(x(i,:));
            x(i,ind)=new;
        end
        xtrain=x(:,1:n_train);
        xtest=x(:,n_train+1:end);
    otherwise
        error('params.scale_meth invalid')
end

% add fixed grid points so we can see shape of estimated set and/or
% estimate its volume
if n_grid > 0
    a0=(1-wf)/2;
    a1=a0+wf;
    h=wf/n_grid;
    a=linspace(a0+h/2,a1-h/2,n_grid);   % 1-d grid
    
    % generate grid
    for k=d-1:-1:0
        xtest_grid(d-k,:)=repmat(reshape(repmat(a,n_grid^k,1),1,n_grid^(k+1)),1,n_grid^(d-k-1));
    end
    xtest=[xtest, xtest_grid];
    if size(xtest,2) ~= n_test + n_grid^d
        warning('grid generation error')
        keyboard
    end
    n_test = size(xtest,2);
    
end

%    x1=reshape(repmat(linspace(0,1,n_grid),[n_grid,1]),1,n_grid^2);
%    x2=repmat(linspace(0,1,n_grid),[1 n_grid]);
%    plot_data=wf*[x1; x2] + (1-wf)/2;
%    xtest=[xtest, plot_data];
%    if size(xtest,2) ~= n_test + n_grid^2
%        warning('grid generation error')
%        keyboard
%    end
%    n_test = size(xtest,2);
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate shifts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(shift_type,'rand')    % random
    sh_list = shift_frac*(2*rand(rf_dim,n_shift)-1);
else    % fixed
    error('fixed shifts not implemented in this version')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% generate "balanced" feature sets for random forests
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_forest > 1
    dim_tot=0;
    feat=[];
    while dim_tot < n_forest*rf_dim
        feat=[feat, randperm(d)];
        dim_tot=dim_tot + d;
    end
    feat=reshape(feat(1:n_forest*rf_dim),[rf_dim,n_forest])';
else
    feat=1:d;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% main loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_vote = n_forest*n_shift;
tr_votes = zeros(n_vote,n_train); % record votes for each test point
ts_votes = zeros(n_vote,n_test); % record votes for each test point

c_lo=0;     % bisection search parameters
c_hi=1;
for outer_iter=1:max_outer
    switch params.prob
        case {'NP','MM'}
            emperr=zeros(1,num_classes);
            bound=zeros(1,num_classes);
        case {'MV','DL'}
            emperr=zeros(1,2);
            bound=0;            % not computed
        otherwise
            emperr=0;
            bound=0;
    end

    ctr=0;  % vote counter
    for rf_ctr=1:n_forest
      if n_forest > 1
          fprintf('random forest counter = %d\n',rf_ctr)
      end
      for sh_ctr=1:n_shift
        if n_shift > 1
          fprintf('  shift counter = %d\n',sh_ctr)
        end
        ctr = ctr+1;
        idx = feat(rf_ctr,:); % features to train on    
        sh = sh_list(:,sh_ctr);
        xtrain_sh = mod(xtrain(idx,:)+repmat(sh,[1,n_train]),1);
        xtest_sh = mod(xtest(idx,:)+repmat(sh,[1,n_test]),1);

        expermode=1-display_dict;

        switch prob
            case {1,2,7}  % PE, CS, RL
                cst_pen = 0;
                [tr_lab,ts_lab,err,pen,sz,dp]=ddt_core(xtrain_sh,ytrain,xtest_sh,max_cuts(idx),1,...
                   'penalty',pen_num,'pen_wt',pen_wt,'costs',costs,'expermode',expermode,'cst_pen',cst_pen);
                %emperr = emperr+err(1);
                %bound = bound + err(1) + pen(1);
                tr_votes(ctr,:)=tr_lab;
                ts_votes(ctr,:)=ts_lab;
            case 3  % NP            
                if shift_inside == 0
                    c_lo=0;     % bisection search parameters
                    c_hi=1;
                end
                for inner_iter=1:max_inner
                    c=(c_lo + c_hi)/2;
                    lam=c/(1-c);

                    % adjust costs and pen weight for NP problem

                    if reweight_meth==2   % treat errors separately
                        costs(i0) = lam/n0;
                        costs(i1) = 1/n1;
                        mod_pen_wt = pen_wt;
                        mod_pen_wt(1)=lam*(1-nu)*pen_wt(1);
                        [tr_lab,ts_lab,err,pen]=ddt_core(xtrain_sh,ytrain,xtest_sh,max_cuts(idx),2,...
                            'penalty',pen_num,'pen_wt',mod_pen_wt,'costs',costs,'expermode',expermode);
                        pf = err(2)/lam;
                        pm = err(3);
                        penf=pen(2);

                    else    % treat errors all at once using cost sensitive classification
                        costs(i1)=1/(n1+lam*n0);
                        costs(i0)=lam/(n1+lam*n0);
                        [tr_lab,ts_lab,err,pen]=ddt_core(xtrain_sh,ytrain,xtest_sh,max_cuts(idx),1,...
                            'penalty',pen_num,'pen_wt',pen_wt,'costs',costs,'expermode',expermode);
                        pf = err(2)*(n1 + lam*n0)/(lam*n0);
                        pm = err(3)*(n1 + lam*n0)/n1;
                        penf=0;
                    end
                    if shift_inside == 0
                        if display_iter == 1
                            fprintf('    iter = %d, lam = %3.4f, pf = %1.4f, pm = %1.4f\n',inner_iter,lam,pf,pm)
                        end

                        if abs(pf-nu*penf-alpha) < tol
                          break
                        end
                        if pf < alpha +nu*penf
                            c_hi = c;
                        else
                            c_lo = c;
                        end
                        expermode=1;
                    end
                end
                tr_votes(ctr,:)=tr_lab;  
                ts_votes(ctr,:)=ts_lab;  
                %emperr = emperr + [pf pm];
            case 4  % MM (two classes only so far)
                warning('minmax algorithm under construction')
                if shift_inside == 0
                    c_lo=0;     % bisection search parameters
                    c_hi=1;
                end
                for inner_iter=1:max_inner
                    c=(c_lo + c_hi)/2;
                    lam=c/(1-c);

                    % adjust costs                
                    if reweight_meth==2   % treat errors separately
                        costs(i0) = lam/n0;
                        costs(i1) = 1/n1;
                        mod_pen_wt = pen_wt;
                        mod_pen_wt(1)=lam*(2-nu)*pen_wt(1);
                        [tr_lab,ts_lab,err,pen]=ddt_core(xtrain_sh,ytrain,xtest_sh,max_cuts(idx),2,...
                            'penalty',pen_num,'pen_wt',mod_pen_wt,'costs',costs,'expermode',expermode);
                        pf = err(2)/lam;
                        pm = err(3);
                        penf = pen(2);
                        penm = pen(3);

                    else    % treat errors all at once using cost sensitive classification
                        costs(i1)=1/(n1+lam*n0);
                        costs(i0)=lam/(n1+lam*n0);
                        [tr_lab,ts_lab,err,pen]=ddt_core(xtrain_sh,ytrain,xtest_sh,max_cuts(idx),1,...
                            'penalty',pen_num,'pen_wt',pen_wt,'costs',costs,'expermode',expermode);
                        pf = err(2)*(n1 + lam*n0)/(lam*n0);
                        pm = err(3)*(n1 + lam*n0)/n1;
                        penf=0; penm=0;
                    end
                    if shift_inside == 0
                        if display_iter == 1
                            fprintf('    iter = %d, lam = %3.4f, pf = %1.4f, pm = %1.4f\n',inner_iter,lam,pf,pm)
                        end

                        if abs(pf+penf-rho*(pm+penm)) < tol
                          break
                        end
                        if pf+penf < rho*(pm+penm)
                            c_hi = c;
                        else
                            c_lo = c;
                        end
                        expermode=1;
                    end
                end
                tr_votes(ctr,:)=tr_lab;  
                ts_votes(ctr,:)=ts_lab;  
                %emperr = emperr + [pf pm];
            case 5  % MV
                if shift_inside == 0
                    c_lo=0;     % bisection search parameters
                    c_hi=1;
                end
                for inner_iter=1:max_inner
                    c=(c_lo + c_hi)/2;
                    lam=c/(1-c);

                    % adjust costs and pen weight
                    costs=(lam/n_train)*ones(1,n_train);
                    %mod_pen_wt=(1-nu*(1-lam))*pen_wt;  % from paper
                    mod_pen_wt=lam*(2-nu)*pen_wt;  % based on NP score type bound
                    [tr_lab,ts_lab,err,pen,sz,dp]=ddt_core(xtrain_sh,ytrain,xtest_sh,max_cuts(idx),3,...
                        'penalty',pen_num,'pen_wt',mod_pen_wt,'costs',costs,'expermode',expermode);
                    mass = 1-err(2)/lam;
                    vol = err(1);

                    if shift_inside == 0
                        if display_iter == 1
                            fprintf('    iter = %d, lam = %3.4f, mass = %1.4f, vol = %1.4f\n',inner_iter,lam,mass,vol)
                        end

                        if abs(mass+nu*pen(2)-alpha) < tol
                          break
                        end
                        if mass >= alpha - nu*pen(2) 
                            c_hi = c;
                        else
                            c_lo = c;
                        end
                        expermode=1;
                    end
                end
                tr_votes(ctr,:)=tr_lab;  
                ts_votes(ctr,:)=ts_lab;  
                %emperr = emperr + [mass vol];
            case 6  % DL
                error('DL Not yet implemented')

        end % of switch statement

      end % loop over shifts
    end % loop over feature sets in random forest
    
    % check stopping criterion
    if shift_inside == 1
        switch prob
            case 3 % NP
                tr_labels = (sum(tr_votes,1)>n_vote/2);
                pf = sum(tr_labels(i0))/n0;
                pm = sum(~tr_labels(i1))/n1;
                if display_iter == 1
                    fprintf('    iter = %d, lam = %3.4f, pf = %1.4f, pm = %1.4f\n',outer_iter,lam,pf,pm)
                end

                if abs(pf-alpha) < tol
                  break
                end
                if pf < alpha 
                    c_hi = c;
                else
                    c_lo = c;
                end
            case 4 % MM: not yet implemented
            case 5 % MV
                mass = sum(sum(tr_votes,1)<=n_vote/2)/n_train;
                if display_iter == 1
                    fprintf('    iter = %d, lam = %3.4f, mass = %1.4f\n',outer_iter,lam,mass)
                end

                if abs(mass-alpha) < tol
                  break
                end
                if mass >= alpha  
                    c_hi = c;
                else
                    c_lo = c;
                end
                expermode=1;
        end
    end
    
end % outer loop over lagrange multiplier
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Count votes
%
% Does not handle multiple classes yet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

labels = (sum(ts_votes,1)>n_vote/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display image of estimated set for 2-d data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if n_grid>0
    set_map=labels(n_test-n_grid^d+1:end);
    if display_est==1
        lims = wf*[0 1]+(1-wf)/2;
        switch d
            case 2
                set_im = reshape(set_map,[n_grid,n_grid]);
                imagesc(lims,lims,.3*set_im+.7, [0,1]);
                colormap(gray)
                axis xy
                hold on
                if prob==5 | prob==6    % give inliers and outliers different labels
                    i0=find(labels(1:n_train)==0); i1=find(labels(1:n_train)==1);
                else % give the different classes different labels
                    i0=find(ytrain==0); i1=find(ytrain==1);
                end

                x_inliers = xtrain(:, i1);
                k = convhull(x_inliers(1, :), x_inliers(2, :));

                scatter(xtrain(1,i0),xtrain(2,i0), 10, 'filled', 'o', 'MarkerFaceColor', [1, 0.6, 0.6]);
                plot(x_inliers(1, k), x_inliers(2, k), 'color', ...
                     [0.4, 0.4, 1], 'LineWidth', 1);
                scatter(xtrain(1,i1),xtrain(2,i1), 10, 'filled', 'o', 'MarkerFaceColor', [0.4, 0.4, 1]);
            case 3
                if prob==5 | prob==6    % give inliers and outliers different markers
                    i0=find(labels(1:n_train)==0); i1=find(labels(1:n_train)==1);
                else % give the different classes different markers
                    i0=find(ytrain==0); i1=find(ytrain==1);
                end
                scatter3(xtrain(1,i0),xtrain(2,i0),xtrain(3,i0),20,'o','filled');
                axis([lims lims lims])
                hold on
                scatter3(xtrain(1,i1),xtrain(2,i1),xtrain(3,i1),100,'rx');
        end
    end
    if display_vol==1
        grid_vol=wf^2*sum(set_map==0)/length(set_map);
    end
    
    % return labels to normal
    labels =labels(1:n_test-n_grid^d);
end

% return statistics
if nargout >= 2
    stats.e_time = etime(clock,t0);
    if prob==5 | prob==6 % MV or DL
        stats.mass = emperr(1)/n_vote;
        if display_vol == 1
            if n_shift == 1
                stats.vol = emperr(2)/n_vote;
            else 
                stats.vol = grid_vol;
            end
        else
            stats.vol = 'NA';
        end
    else
        stats.emperr = emperr/n_vote;
    end
    stats.bound = bound/n_vote;
    stats.sh = sh_list;
    stats.feat = feat;
    stats.depth = [];
    stats.size = [];
end