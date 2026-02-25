function [ sig_denoised ] = cadzow_filter( input_args, r)
% CADZOW_FILTER Summary of this function goes here
% Detailed explanation goes here
% both input and output are one dimensional signals
% Author: Xi Peng
H = form_Hankel( input_args );
[u s v]=svd(H, 'econ');

s1=zeros(size(s));
for i=1:r
    s1(i,i)=s(i,i);
end
H_denoised=u*s1*v';
sig_denoised = extract_1Dsignal_from_Hankel( H_denoised );% average 
end

