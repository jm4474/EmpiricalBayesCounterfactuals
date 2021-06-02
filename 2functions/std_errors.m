function [standard_errors,Asy_Cov] = std_errors(Y,X,Z,...
                                              beta,gamma,...
                                              hessian_beta,...
                                              hessian_gamma)
%% std_errors: This function calculates standard errors for parameters
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [standard_errors,Asy_Cov] = std_errors(Y,X,Z,...
%                                              beta,gamma,...
%                                              hessian_beta,...
%                                              hessian_gamma)
%
% INPUT:
% Y: matrix of size J x T
% X: Array of size  J x T x K 
% Z: Array of size  J x L
% beta: K x 1 parameter vector
% gamma: (L+1) x 1 parameter vector
% hessian_beta: The hessian matrix of the objective function used to
% estimate beta
% hessian_gamma: The hessian matrix of the objective function used to
% estimate gamma

%% Dimensions

J = size(Y,1);

T = size(Y,2);

K = size(X,3);

L = size(Z,2); 

Zaux ...
  = [ones(J,1),Z]; 

%% Sandwich covariance matrix


bread ...
    = zeros(K+L+1,K+L+1,size(Y,1));

two_step_correction ...
    = zeros(K,L+1,size(Y,1));

for j_LEA = 1: size(Y,1)
   
    %Derivative of first moment function w.r.t beta
    
    Xaux_j = reshape(X(j_LEA,:,:),[T,K]);
             %[X_{j1},...,X_{jK}], where X_{jk} is the time series of cov k
    
    aux1   = -sum(Xaux_j'.*(Y(j_LEA,:)),2)./T; 
             %K x 1 vector equal to (1/T) sum_{t=1}^{T} X_{jt} Y_{jt}
    
    aux2   = exp(Xaux_j*beta);
             %T x 1 vector containing exp(X_{jt}'beta)
    
    aux3   = sum(Y(j_LEA),2).*aux2./sum(aux2,1);
             %T x 1 vector containing bar(Y_j)*exp(X_{jt}'beta)/sum_{t=1}^Texp(X_{jt}'beta) 
    
    DM1Dbeta   ...
           = aux1 +  (sum(Xaux_j.*aux3,1)'./T);
             %K x 1 vector containing the derivative of the 
             %first moment condition evaluated at beta
    
    %Derivative of second moment function w.r.t beta
    
    Zaux_j = Zaux(j_LEA,:);
    
    aux4   = sum(Y(j_LEA),2).*aux2./(sum(aux2,1).^2);
             %T x 1 vector containing
             %bar(Y_j)*exp(X_{jt}'beta)/(sum_{t=1}^Texp(X_{jt}'beta))^2
    
    two_step_correction(:,:,j_LEA) ...
           = (sum(Xaux_j.*aux4,1)').*exp(Zaux_j*gamma)*Zaux_j;
             %K times (L+1)
             
    %Derivative of second moment function w.r.t gamma
    
    DM2Dgamma ...
           = -(sum(Y(j_LEA),2)./(sum(aux2,1))...
           - exp(Zaux_j*gamma)).*exp(Zaux_j*gamma)*Zaux_j';
           % L+1 x 1 vector of derivatives
    
    bread(:,:,j_LEA)...
           = [DM1Dbeta;DM2Dgamma]*[DM1Dbeta;DM2Dgamma]';
       
    clear Xauxj Zaux_j aux1 aux2 aux3 aux4 DM1Dbeta DM2Dgamma
             
end

%% 

hessian_all ...
    = [hessian_beta , zeros(K,L+1) ; ...
       mean(two_step_correction,3)', hessian_gamma];

Asy_Cov ...
    = (full(hessian_all)^(-1))*mean(bread,3)*(full(hessian_all)'^(-1));

J   = size(Y,1);

standard_errors ...
    = diag(Asy_Cov./J).^.5;

end

