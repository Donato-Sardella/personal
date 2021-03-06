function [J, grad] = cofiCostFunc(params, Y, R, num_users, num_movies, ...
                                  num_features, lambda)
%COFICOSTFUNC Collaborative filtering cost function
%   [J, grad] = COFICOSTFUNC(params, Y, R, num_users, num_movies, ...
%   num_features, lambda) returns the cost and gradient for the
%   collaborative filtering problem.
%

% Unfold the U and W matrices from params
X = reshape(params(1:num_movies*num_features), num_movies, num_features);
Theta = reshape(params(num_movies*num_features+1:end), ...
                num_users, num_features);

            
% You need to return the following values correctly
J = 0;
X_grad = zeros(size(X));
Theta_grad = zeros(size(Theta));

% ====================== YOUR CODE HERE ======================
% Instructions: Compute the cost function and gradient for collaborative
%               filtering. Concretely, you should first implement the cost
%               function (without regularization) and make sure it is
%               matches our costs. After that, you should implement the 
%               gradient and use the checkCostFunction routine to check
%               that the gradient is correct. Finally, you should implement
%               regularization.
%
% Notes: X - num_movies  x num_features matrix of movie features
%        Theta - num_users  x num_features matrix of user features
%        Y - num_movies x num_users matrix of user ratings of movies
%        R - num_movies x num_users matrix, where R(i, j) = 1 if the 
%            i-th movie was rated by the j-th user
%
% You should set the following variables correctly:
%
%        X_grad - num_movies x num_features matrix, containing the 
%                 partial derivatives w.r.t. to each element of X
%        Theta_grad - num_users x num_features matrix, containing the 
%                     partial derivatives w.r.t. to each element of Theta
%


%Y: a matrix of movie ratings: Dimensions are (movies 1682 x users 944)

%X: a matrix of movie features (0 to 5): Dimensions are (movies 1682 x features 10) 

%Theta: a matrix of feature weights: Dimensions are (users 944 x features 10)

%Error factor= X*Theta'=(1682*10) X (10*944)=(1682*944)

J=0.5*sum(sum((R.*((X*Theta')-Y).^2)));
reg_J=((lambda/2)*(sum(sum(Theta.^2))))+((lambda/2)*(sum(sum(X.^2))));
J=J+reg_J;
%The X gradient is the product of the error factor and the Theta matrix. 
%The sum is computed automatically by the vector multiplication. Dimensions are (movies 1682 x features 10) 

X_grad=(R.*((X*Theta')-Y))*Theta;
X_grad=X_grad+(lambda*X);


%The Theta gradient is the product of the error factor and the X matrix. A transposition may be needed. 
%The sum is computed automatically by the vector multiplication. Dimensions are (users 944 x features 10) 

Theta_grad=(R.*((X*Theta')-Y))'*X;
Theta_grad=Theta_grad+(lambda*Theta);

% =============================================================

grad = [X_grad(:); Theta_grad(:)];

end
