% This is an implementation of the robust PCA method as described in 
% 
% A. Podosinnikova, S. Setzer, and M. Hein,
% Robust PCA: Optimization of the Robust Reconstruction Error 
% over the Stiefel Manifold, GCPR 2014
% 
% All copyrights remain with the authors - you have to cite the above
% paper when using this code.
% 
%
% With any problems/questions/suggestions/success stories contact:
%
% Matthias Hein - hein@cs.uni-sb.de
% or
% Anastasia Podosinnikova - anastasia.podosinnikova@inria.fr
% 
% 
% Input:
% X is a n by p data matrix: 
%       n - the number of observations (in rows)
%       p - the dimension of each observation
% param is a struct of parameters:
%       param.k - the number of robust PCs to recover
%       param.t - a lower bound on the number of true observations (def. n/2)
%       param.nruns - the number of random restarts .............. (def. 10)
%       param.eps - solution tolerance ........................... (def. 1E-5)
%       param.maxiter - maximal number of iterations ............. (def. 1E5)
%       param.isdisplay - display objectives on the screen if 1 .. (def. 0)
%
%
% Default usage:
%           [U,m] = trpca(X,k);
%           where k is the number of desired robust PCs
%
% Otherwise, e.g.
%           param.k = k; param.nruns = nruns; param.t = t;
%           [U,m] = trpca(X,param);
%
% Output:
%           U is a p by k matrix of the robust PCs (PCs are in columns)
%           m is a p column vector of the robust center
%
% Advice on the choice of nruns, i.e. the number of random restarts. Turn
% on isdisplay option to track the objective values achieved after every
% random restart and, if there are a lot of different values and no
% repetitions (hence, too many local extreema), increase the nruns value.
%
function [U,m] = trpca(X,param)
  [p,d,t,nruns,maxiter,eps,isdisplay] = initialize_parameters(X,param);
  objs=zeros(nruns,1); Us=cell(nruns,1); ms=cell(nruns,1);
  for i=1:nruns
    U0=sample_from_Stiefel_manifold(p,d);
    [U,m,obj] = trpca_one_run(X,U0,t,maxiter,eps);
    objs(i)=obj; Us{i}=U; ms{i}=m;
    if isdisplay, fprintf('run = %3i \t objective obj = %3.9f\n',i,obj); end
  end
  [~,imin]=min(objs); U=Us{imin}; m=ms{imin};
end


function [U,m,obj] = trpca_one_run(X,U,t,maxiter,eps)
  n = size(X,1);
  m = median(X);
  objold = 2^53-1; % bitmax
  iter=1;
  while 1>0
    % (1) update robust PCs U
    XC = X - ones(n,1)*m;
    normXCs = sum(XC.^2,2);
    proj = (XC*U); 
    r = normXCs - sum(proj.*proj,2); 
    [sr,ir] = sort(r,'ascend');
    obj = sum(sr(1:t));
    
    U = update_U(XC,t,U,ir,obj,normXCs,maxiter,eps);
    
    % (2) update robust center m
    proj = (XC*U); 
    r = normXCs - sum(proj.^2,2); 
    [sr,ir] = sort(r,'ascend');
    m = 1/t * sum(X(ir(1:t),:));
    obj = sum(sr(1:t));
    
    % stopping criterion
    if objold - obj <= eps, break, end
    objold = obj;
    if iter >= maxiter, warning('Maxiter:main','Maxiter:main'), break, end
    iter=iter+1;
  end
end


function [U,obj] = update_U(XC,t,U,ir,obj,normXC2,maxiter,eps)
  oldobj = 2^53-1; % bitmax
  iter=1; 
  while 1>0
    % compute gradient
    x=XC(ir(1:t),:); grad=x'*(x*U);
    
    % compute polar decomposition
    [W,~,V]=svd(grad,0);
    U=W*V';
    
    % compute objective
    proj=XC*U;
    r = normXC2 - sum(proj.^2,2);
    [sr,ir]=sort(r,'ascend');
    obj=sum(sr(1:t));
    
    % stopping criterion
    if oldobj - obj <= eps, break, end
    oldobj = obj;
    if iter >= maxiter, warning('Maxiter:U','Maxiter:U'); break, end
    iter=iter+1;
  end
end


function U = sample_from_Stiefel_manifold(p,k)
  X = randn(p,k); U = orth(X);
end


function [p,k,t,nruns,maxiter,eps,isdisplay] = initialize_parameters(X,param)
  % n - number of observations in X
  % p - dimension (number of features) of an observation in X
  [n,p]=size(X);
  if ~isstruct(param), 
    k = param; clear param;
    param.k = k;
  end
  
  % k - number of PCs to recover
  if isfield(param,'k')
    k = param.k;
    if (k<1) || (k>p) || (floor(k)~=k), error('Bad value of k'), end
  else
    k = min(10,ceil(p/2));
    warning('no:k','Number of robust PCs is not provided');
  end;
  
  % t - lower bound on the number of true observations (
  if isfield(param,'t')
    t = param.t;
    if (t<ceil(n/2)) || (t>n) || (floor(t)~=t), error('Bad value of t'), end
  else
    t = ceil(n/2);
  end
  
  % nruns - number of random restarts of TRPCA
  if isfield(param,'nruns')
    nruns = param.nruns;
    if (nruns<0) || (floor(nruns)~=nruns), error('Bad value of nruns'), end
  else
    nruns = 10;
  end
  
  % eps - tolerance
  if isfield(param,'eps')
    eps = param.eps;
    if eps<0, error('Bad value of eps'), end
  else
    eps = 1E-5;
  end
  
  % maxiter - maximal number of iterations
  if isfield(param,'maxiter')
    maxiter = param.maxiter;
    if maxiter < 0, error('Bad value of maxiter'), end
  else
    maxiter = 1E5;
  end
  
  % isdisplay - display objectives on the screen if 1
  if isfield(param,'isdisplay')
    isdisplay = param.isdisplay;
    if isdisplay~=1, isdisplay = 0; end
  else
    isdisplay = 0;
  end
  
end


