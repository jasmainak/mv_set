%close all;
h = figure('Name', label, 'NumberTitle', 'off');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set parameters

params.max_cuts=8;   
params.penalty = penalty;
params.pen_wt = 1;
params.n_shift = 11;   
params.shift_frac = .1;  
params.display_est = 1;
params.display_vol = 1;
params.display_iter = 1;
params.shift_inside = 1;
params.width_frac = .8;
params.reweight_meth = 1;
params.max_iter = 15;
%params.n_grid = 400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in data. Change the switch argument to select a different dataset.
switch 5
    case 1 
        %params.prob = 'PE';
        %params.prob = 'NP';
        params.prob = 'MV';    % be sure to reset params.alpha
        
        params.display_est=1;
        
        load banana 
        
        if strcmp(params.prob,'NP')
            params.alpha=0.9;
            %xtest=xtrain;
            %ytest=ytrain;
        end
        if strcmp(params.prob,'MV') 
            params.alpha=0.99;
            xtrain=xtrain(:,find(ytrain==0));
            xtest=xtest(:,find(ytest==0));
            n_train=size(xtrain,2);
            n_test=size(xtest,2);
            ytrain=zeros(1,size(xtrain,2));
        end
        
    case 2  % gaussian bumps with spherical variance and different means
        % tweak p below to change ratio of positives to negatives
        
		%params.prob = 'PE';  % PE favors majority class 
        params.prob = 'NP';   % but NP should not reflect class freqs.

        n_train=500;
        n_test=1000;
        d=2; % dimension
        k=1; % difference in means
        p=.3; % proportion of class 0
        mu0=zeros(d,1);  mu1=k*ones(d,1);
        cov0=eye(d);  cov1=eye(d);
        % generate train data
        n0=floor(n_train*p); n1=n_train - n0;
        x0=mvg(n0,d,mu0,cov0);  x1=mvg(n1,d,mu1,cov1);
        xtrain=[x0,x1];
        ytrain=ones(1,n_train);        ytrain(1:n0)=0;
        % generate test data
        n0=floor(n_test*.5);
        n1=n_test - n0;
        x0=mvg(n0,d,mu0,cov0);        x1=mvg(n1,d,mu1,cov1);
        xtest=[x0,x1];
        ytest=ones(1,n_test);        ytest(1:n0)=0;

        if strcmp(params.prob,'NP') % determine optimal pm
            params.alpha=0.1;
            snr = d*k^2;    % variance = 1
            pm_opt = .5+.5*erf((sqrt(2)*erfinv(1-2*params.alpha)-sqrt(snr))/sqrt(2));
            fprintf('optimal miss probability is %1.4f\n',pm_opt);
        end

    case 3 % spherical gaussian
        params.prob = 'MV';  
        params.alpha = 0.9;
        
        n_train=500;
        n_test=1000;
        d=3;
        mu=zeros(d,1);
        covarm=eye(d);
        xtrain=mvg(n_train,d,mu,covarm);
        xtest = xtrain; n_test = n_train;
%        xtest=mvg(n_test,d,mu,covarm);
        ytrain=zeros(1,size(xtrain,2));
    case 4 % 2d mixture data used in paper
        params.prob = 'MV';  
        params.alpha = 0.9;

        n_train=500;
        n_test=500;
        xtrain = [mvg(n_train/2,2,[0 0]', [1 .7; .7 1]), ...
            mvg(n_train/2,2,[0 -2]', [1 -.95; -.95 1])];
        xtest = [mvg(n_test/2,2,[0 0]', [1 .7; .7 1]), ...
            mvg(n_test/2,2,[0 -2]', [1 -.95; -.95 1])];
        ytrain=zeros(1,size(xtrain,2));
        xtest = xtrain;
    case 5 % new
        params.prob = 'MV';
        params.alpha = 0.9;
        
        params.display_est=1;
       
        fid = fopen(fname);
        C = textscan(fid, '%s%f%f');
        x = [C{2}, C{3}];
        last_train = floor(size(x, 1) / 2);
        xtrain = x(1:last_train, :)'; xtest = x(last_train:end, :)';
        ytrain = zeros(1, size(xtrain, 2));
        
        n_train = size(xtrain, 2); n_test = size(xtest, 2);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Learn set estimate and predict labels

[labels,stats] = ddt_learn(xtrain,ytrain,xtest,params);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Report results

prob=params.prob;
switch prob
    case {'PE','CS','NP','RL','MM'}
        if prob=='RL'
           ytest = (ytest > params.level); 
        end
        pe=sum(labels~=ytest)/length(ytest);
        fpos=sum((labels~=ytest) .* (ytest==0))/sum(ytest==0);
        fneg=sum((labels~=ytest) .* (ytest==1))/sum(ytest==1);

        disp(['Misclassification rate = ' num2str(pe)])
        disp(['False positive rate = ' num2str(fpos)])
        disp(['False negative rate = ' num2str(fneg)])
    case {'MV','DL'}
        % depends on what you want to do
        emass=sum(labels(1:n_test)==0)/n_test;
        disp(['Test mass = ' num2str(emass)])
        if ~strcmp(stats.vol,'NA')
            vol=stats.vol;
            disp(['Volume = ' num2str(vol)])
        end

end