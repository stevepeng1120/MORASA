function [ output_args ] = form_Hankel( input_args )
% FORM_HANKEL Summary of this function goes here
%   Detailed explanation goes here
% input_args is the linear predictable signal
% Author: Xi Peng
sig = reshape(input_args,1,[]);
N = length(sig);
Ny = ceil(N/2);
Nx = N-Ny+1;
H = zeros(Nx,Ny);
for i = 1:Nx
    H(i,:) = sig(i:i+Ny-1);
end
output_args = H;
end

