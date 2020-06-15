% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------------------------------------
% Created by isardSAT S.L. 
% --------------------------------------------------------
%
% CryoSat 2 calibration over transponders
% 
% This code implements the algorithm as described in the
% ISARD_ESA_CR2_TRP_CAL_DPM_030 2.b of 26/05/2011
%
% ---------------------------------------------------------
% FMINSEARCH_STACKCAK: code based on the Matlab FMINSEARCH code.
% The number of iterations and the tolerances have been increased.
% Some parameter values have been modified
% ----------------------------------------------------------
% 
% Author:   Mercedes Reche / Pildo Labs
% Reviewer: Monica Roca / isardSAT
%
%
% $Id: fminsearch_stackcal.m 137 2005-12-20 07:31:57Z merche $
%
% This software is subject to the conditions 
% set forth in the ESA contract CCN2 of contract #C22114 / #4200022114
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,fval,exitflag,output] = fminsearch_stackcal(funfcn,x,options,varargin)
%FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%   X = FMINSEARCH(FUN,X0) starts at X0 and finds a local minimizer X of the
%   function FUN. FUN accepts input X and returns a scalar function value 
%   F evaluated at X. X0 can be a scalar, vector or matrix. 
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, and MaxIter. 
%
%   X = FMINSEARCH(FUN,X0,OPTIONS,P1,P2,...) provides for additional
%   arguments which are passed to the objective function, F=feval(FUN,X,P1,P2,...).
%   Pass an empty matrix for OPTIONS to use the default values.
%   (Use OPTIONS = [] as a place holder if no options are set.)
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns a string EXITFLAG that 
%   describes the exit condition of FMINSEARCH.  
%   If EXITFLAG is:
%     1 then FMINSEARCH converged with a solution X.
%     0 then the maximum number of iterations was reached.
%   
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value 
%     SIN evaluated at X.
%
%     FUN can also be an inline object:
%        X = fminsearch(inline('norm(x)'),[1;2;3])
%     returns a minimum near [0;0;0].
%
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, @, INLINE.

%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1): 
%   p.112-147, 1998.

%   Copyright 1984-2002 The MathWorks, Inc. 
%   $Revision: 1.1 $  $Date: 2005-12-20 08:31:57 +0100 (mar, 20 dic 2005) $


defaultopt = struct('Display','notify','MaxIter','2000*numberOfVariables',...
   'MaxFunEvals','2000*numberOfVariables','TolX',1e-12,'TolFun',1e-12);

% If just 'defaults' passed in, return the default options in X
if nargin==1 && nargout <= 1 && isequal(funfcn,'defaults')
   x = defaultopt;
   return
end

if nargin < 2, 
   error('FMINSEARCH requires at least two input arguments'); 
end

if nargin<3, options = []; end
n = prod(size(x));
numberOfVariables = n;

printtype = optimget(options,'Display',defaultopt,'fast');
tolx = optimget(options,'TolX',defaultopt,'fast');
tolf = optimget(options,'TolFun',defaultopt,'fast');
maxfun = optimget(options,'MaxFunEvals',defaultopt,'fast');
maxiter = optimget(options,'MaxIter',defaultopt,'fast');

% In case the defaults were gathered from calling: optimset('fminsearch'):
if ischar(maxfun)
   if isequal(lower(maxfun),'2000*numberofvariables')
      maxfun = 2000*numberOfVariables;
   else
      error('Option ''MaxFunEvals'' must be an integer value if not the default.')
   end
end
if ischar(maxiter)
   if isequal(lower(maxiter),'2000*numberofvariables')
      maxiter = 2000*numberOfVariables;
   else
      error('Option ''MaxIter'' must be an integer value if not the default.')
   end
end

switch printtype
case 'notify'
   prnt = 1;
case {'none','off'}
   prnt = 0;
case 'iter'
   prnt = 3;
case 'final'
   prnt = 2;
case 'simplex'
   prnt = 4;
otherwise
   prnt = 1;
end

header = ' Iteration   Func-count     min f(x)         Procedure';

% Convert to inline function as needed.
funfcn = fcnchk(funfcn,length(varargin));

n = prod(size(x));

% Initialize parameters
rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
onesn = ones(1,n);
two2np1 = 2:n+1;
one2n = 1:n;

% Set up a simplex near the initial guess.
xin = x(:); % Force xin to be a column vector
v = zeros(n,n+1); fv = zeros(1,n+1);
v(:,1) = xin;    % Place input guess in the simplex! (credit L.Pfeffer at Stanford)
x(:) = xin;    % Change x to the form expected by funfcn 
fv(:,1) = feval(funfcn,x,varargin{:}); 

% Following improvement suggested by L.Pfeffer at Stanford
usual_delta = 0.05;           % 5 percent deltas for non-zero terms
zero_term_delta = 0.00025;   % Even smaller delta for zero elements of x
for j = 1:n
   y = xin;
   if y(j) ~= 0
      y(j) = (1 + usual_delta)*y(j);
   else 
      y(j) = zero_term_delta;
   end  
   v(:,j+1) = y;
   x(:) = y; f = feval(funfcn,x,varargin{:});
   fv(1,j+1) = f;
end     

% sort so v(1,:) has the lowest function value 
[fv,j] = sort(fv);
v = v(:,j);

how = 'initial';
itercount = 1;
func_evals = n+1;
if prnt == 3
   disp(' ')
   disp(header)
   disp([sprintf(' %5.0f        %5.0f     %12.6g         ', itercount, func_evals, fv(1)), how]) 
elseif prnt == 4
   clc
   formatsave = get(0,{'format','formatspacing'});
   format compact
   format short e
   disp(' ')
   disp(how)
   v
   fv
   func_evals
end
exitflag = 1;

% Main algorithm
% Iterate until the diameter of the simplex is less than tolx
%   AND the function values differ from the min by less than tolf,
%   or the max function evaluations are exceeded. (Cannot use OR instead of AND.)
while func_evals < maxfun && itercount < maxiter
%while func_evals < 300000 & itercount < maxiter%A�ADIDO PARA PROBAR
    if max(max(abs(v(:,two2np1)-v(:,onesn)))) <= tolx && ...
         max(abs(fv(1)-fv(two2np1))) <= tolf
      break
   end
   how = '';
   
   % Compute the reflection point
   
   % xbar = average of the n (NOT n+1) best points
   xbar = sum(v(:,one2n), 2)/n;
   xr = (1 + rho)*xbar - rho*v(:,end);
   x(:) = xr; fxr = feval(funfcn,x,varargin{:});
   func_evals = func_evals+1;
   
%    x_save(func_evals)=x;
%    f_save(func_evals)=fxr;
%    
   if fxr < fv(:,1)
      % Calculate the expansion point
      xe = (1 + rho*chi)*xbar - rho*chi*v(:,end);
      x(:) = xe; fxe = feval(funfcn,x,varargin{:});
      func_evals = func_evals+1;
%       x_save(func_evals)=x;
%       f_save(func_evals)=fxe;
      if fxe < fxr
         v(:,end) = xe;
         fv(:,end) = fxe;
         how = 'expand';
      else
         v(:,end) = xr; 
         fv(:,end) = fxr;
         how = 'reflect';
      end
   else % fv(:,1) <= fxr
      if fxr < fv(:,n)
         v(:,end) = xr; 
         fv(:,end) = fxr;
         how = 'reflect';
      else % fxr >= fv(:,n) 
         % Perform contraction
         if fxr < fv(:,end)
            % Perform an outside contraction
            xc = (1 + psi*rho)*xbar - psi*rho*v(:,end);
            x(:) = xc; fxc = feval(funfcn,x,varargin{:});
            func_evals = func_evals+1;
%             x_save(func_evals)=x;
%             f_save(func_evals)=fxc;
            
            if fxc <= fxr
               v(:,end) = xc; 
               fv(:,end) = fxc;
               how = 'contract outside';
            else
               % perform a shrink
               how = 'shrink'; 
            end
         else
            % Perform an inside contraction
            xcc = (1-psi)*xbar + psi*v(:,end);
            x(:) = xcc; fxcc = feval(funfcn,x,varargin{:});
            func_evals = func_evals+1;
%             x_save(func_evals)=x;
%             f_save(func_evals)=fxcc;
            if fxcc < fv(:,end)
               v(:,end) = xcc;
               fv(:,end) = fxcc;
               how = 'contract inside';
            else
               % perform a shrink
               how = 'shrink';
            end
         end
         if strcmp(how,'shrink')
            for j=two2np1
               v(:,j)=v(:,1)+sigma*(v(:,j) - v(:,1));
               x(:) = v(:,j); fv(:,j) = feval(funfcn,x,varargin{:});
            end
            func_evals = func_evals + n;
         end
      end
   end
   [fv,j] = sort(fv);
   v = v(:,j);
   itercount = itercount + 1;
   if prnt == 3
   disp([sprintf(' %5.0f        %5.0f     %12.6g         ', itercount, func_evals, fv(1)), how]) 
   elseif prnt == 4
      disp(' ')
      disp(how)
      v
      fv
      func_evals
   end  
end   % while


x(:) = v(:,1);
if prnt == 4,
   % reset format
   set(0,{'format','formatspacing'},formatsave);
end
output.iterations = itercount;
output.funcCount = func_evals;
output.algorithm = 'Nelder-Mead simplex direct search';

fval = min(fv); 
if func_evals >= maxfun 
   if prnt > 0
%       disp(' ')
%       disp('Exiting: Maximum number of function evaluations has been exceeded')
%       disp('         - increase MaxFunEvals option.')
%       msg = sprintf('         Current function value: %f \n', fval);
%       disp(msg)
   end
   exitflag = 0;
elseif itercount >= maxiter 
   if prnt > 0
      disp(' ')
      disp('Exiting: Maximum number of iterations has been exceeded')
      disp('         - increase MaxIter option.')
      msg = sprintf('         Current function value: %f \n', fval);
      disp(msg)
   end
   exitflag = 0; 
else
   if prnt > 1
      convmsg1 = sprintf([ ...
         '\nOptimization terminated successfully:\n',...
         ' the current x satisfies the termination criteria using OPTIONS.TolX of %e \n',...
         ' and F(X) satisfies the convergence criteria using OPTIONS.TolFun of %e \n'
          ],tolx, tolf);
      disp(convmsg1)
   end
   exitflag = 1;
end




