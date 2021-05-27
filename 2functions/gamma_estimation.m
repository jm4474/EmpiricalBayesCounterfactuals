function [f,g] = gamma_estimation(Y,X,Z,beta,gamma)
%% gamma_estimation: This function estimates 
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [f,g] = gamma_estimation(Y,X,Z,beta,gamma)
% 
% INPUT:
% Y: matrix of size J x T
% X: Array of size  J x T x K 
% Z: Array of size  J x L
% beta:    K x 1 parameter vector
% gamma1: (L+1) x 1 parameter vector

%%
J = size(Y,1);

T = size(Y,2);

K = size(X,3);

L = size(Z,2); 

Zaux ...
  = [ones(J,1),Z];  

%% Create auxiliary outcome variable for nonlinear least squares

aux1 = zeros(J,T);

for t_period = 1:T

aux1(:,t_period) ...
  = exp(reshape(X(:,t_period,:),[J,K])*beta);  
  %J x 1 vector with rows exp([X'_{jt}beta]) 
  
end
 
Y_new ...
  = sum(Y,2)./sum(aux1,2); 
    %J x 1 matrix: ( sum_{t=1}^{T} Y_{jt} ) /  \sum_{t=1}^{T} exp(X_jt'beta)  

%% Objective Function

f ...
  = .5*sum((Y_new - exp( Zaux*gamma )).^2,1)...
    ./(J);

%% Gradient and Hessian

if nargout > 1
    
g ...
  = -(...
     sum(((Y_new -  exp( Zaux*gamma )).*exp(Zaux*gamma)).* Zaux,1)...
     )'...
    ./(J);   

  %if nargout > 2
         
      
      
  %end

end

end

