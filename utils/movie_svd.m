function [U_k, s_k, V_k, denoised_movie] = movie_svd(movie, k)
% denoise a movie with SVD using using k singular values
% default value for k is 500
%
% Usage:
%       denoised_movie = movie_svd(movie, k);

if nargin < 2 || isempty(k), k = 500; end

[X, Y, T] = size(movie);
if T>1
    movie = reshape(double(movie), [X.*Y, T]);
end

[U, s, V] = svd(movie);
U_k = U(:,1:k);
s_k = s(1:k, 1:k);
V_k = V(:,1:k);

denoised_movie = U_k * s_k * V_k';
if T > 1
    denoised_movie = reshape(denoised_movie, [X, Y, T]);
end