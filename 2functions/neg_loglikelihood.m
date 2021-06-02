function [f,g,H] = neg_loglikelihood(Y,X,beta)
%% neg_loglikelihood: This function returns the negative log-likelihood of the multinomial, used in estimating beta. 
% 
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [f,g,H] = neg_loglikelihood(Y,X,beta)
%
% INPUT:
% Y: matrix of size J x T
% X: Array of size  J x T x K 
% beta: K x 1 parameter vector

%% 1) Objective Function

J = size(Y,1);

T = size(Y,2);

K = size(X,3);

aux1 ...
  = zeros(J,K,T); 

aux2 ...
  = zeros(J,T);  

for t_period = 1:T

aux1(:,:,t_period) ...
  = Y(:,t_period).*reshape(X(:,t_period,:),[J,K]); 
  %J x K matrix with rows Y_{jt}X_{jt})'

aux2(:,t_period) ...
  = exp(reshape(X(:,t_period,:),[J,K])*beta);  
  %J x 1 vector with rows exp([X'_{jt}beta]) 
  
end

aux3 ...
  = sum(aux1,3);      %J x K  matrix: \sum_{t=1}^{T} Y_{jt}X_{jt})'

aux4 ...
  = log(sum(aux2,2)); %J x 1 matrix: log(\sum_{t=1}^{T} exp(X_jt'beta)) 

f ...
  = (-(sum(aux3,1)*beta) ...
    + sum(aux4.*sum(Y,2),1))...
    ./(T*J);
% -(sum_{j=1}^{J}\sum_{t=1}^{T} Y_{jt}X_{jt})')*beta  
% + sum_{j=1}^{J} log(\sum_{t=1}^{T} exp(X_jt'beta)) bar(Y)_j
% divided by T*J
%% 2) Gradient and Hessian

if nargout > 1

aux6 ...
  = zeros(J,T,K); 
  
for t_period = 1:T
    
    aux6(:,t_period,:) ...
             = aux2(:,t_period).*X(:,t_period,:);
               %exp(X_jt'beta)*X_{jt}' 
end

aux7 ...
  = reshape(sum(aux6,2),[J,K]); %sum_{t=1}^{T} exp(X_jt'beta)X'_{tj}
    
g ...
  = ( (-sum(aux3,1)') ...
    + sum(((aux7.*sum(Y,2))./sum(aux2,2)),1)')...
    ./(T*J);   

% -(sum_{j=1}^{J}\sum_{t=1}^{T} Y_{jt}X_{jt})) 
% + sum_{j=1}^{J} (sum_{t=1}^{T} exp(X_jt'beta)X'_{tj}*bar(Y)_j / sum_{t=1}^{T}exp([X'_{jt}beta]) ) 

  if nargout > 2
      
      aux8 = zeros(K,K,T);
      
      aux9 = zeros(J,T,K);
      
      aux10 = zeros(K,K,T);
      
      for t_period = 1:T
      
      aux8 (:,:,t_period) = (1/(T*J)).*...
                              reshape(X(:,t_period,:).*...
                                      sum(Y,2)...
                                      ,[J,K])'*...
                              reshape(aux6(:,t_period,:)./...
                                      sum(aux2,2)...
                                      ,[J,K]);
                                  
      aux9(:,t_period,:) ...
                           = (aux2(:,t_period).^2).*X(:,t_period,:);
                           %exp(2X_jt'beta)*X_{jt}' 
                          
                           
      aux10(:,:,t_period) ... 
                           = (1/(T*J)).*...
                              reshape(X(:,t_period,:).*...
                                      sum(Y,2)...
                                      ,[J,K])'*...
                              reshape(aux9(:,t_period,:)./...
                                      (sum(aux2,2).^2)...
                                      ,[J,K]);      
      
      end    
          
      H = (sum(aux8,3) - sum(aux10,3))./(T*J);
      
  end

end

end