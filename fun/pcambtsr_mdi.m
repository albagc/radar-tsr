function [X,m,S,It,diff,Xrec]=pcambtsr_mdi(X,A,maxiter,tol)
% Missing Data Imputation Toolbox v1.0
% A. Folch-Fortuny, F. Arteaga and A. Ferrer
% Copyright (C) 2015 A. Folch-Fortuny and F. Arteaga
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%
% INPUTS
%
% X: data matrix with NaNs for the missing data.
% A: number of principal components.
% maxiter: maximum iterations (e.g. 5000).
% tol: tolerance (e.g. 10e-10).
%
% OUTPUTS
%
% X: original data set with the imputed values.
% m: estimated mean vector of X.
% S: estimated covariance matrix of X.
% It: number of iterations.
% Xrec: PCA reconstruction of X with A components. 


[n,p]=size(X);
for i=n:-1:1
  r=~isnan(X(i,:));
  pat(i).O=find(r==1); % observed variables
  pat(i).M=find(r==0); % missing variables
  pat(i).nO=size(pat(i).O,2); % number of observed variables
  pat(i).nM=size(pat(i).M,2); % number of missing variables
end
mis=isnan(X);
[r, c]=find(isnan(X));
X(mis)=0;
meanc=sum(X)./(n-sum(mis));
for k=1:length(r)
  X(r(k),c(k))=meanc(c(k));
end

diff=100;
It=0;
while and(It<maxiter, diff>tol)
  It=It+1;
  Xmis=X(mis);
  mX=mean(X);
  S=cov(X);
  Xc=X-ones(n,1)*mX;
  if n>p
      [~, ~, V]=svd(Xc,0); 
  else 
      [V, ~, ~]=svd(Xc',0); 
  end
  V=V(:,1:A);
  for i=1:n           % for each row
    if pat(i).nM>0   % if there are missing values
      L=V(pat(i).O,1:min(A,pat(i).nO)); % L is the key matrix
      S11=S(pat(i).O,pat(i).O);
      S21=S(pat(i).M,pat(i).O);
      z1=Xc(i,pat(i).O)';
      z2=S21*L*pinv(L'*S11*L)*L'*z1;
      Xc(i,pat(i).M)=z2';
    end
  end
  X=Xc+ones(n,1)*mX;
  d=(X(mis)-Xmis).^2;
  diff=mean(d);
  
end
S=cov(X);
m=mean(X);
[~, ~, v]=svd(S,0);
P=v(:,1:A);
T=(X-ones(n,1)*m)*P;
Xrec=ones(n,1)*m+T*P';

