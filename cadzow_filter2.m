function [ sig_denoised ] = cadzow_filter2( input_args, r)
%CADZOW_FILTER2 此处显示有关此函数的摘要
%   此处显示详细说明
%CADZOW_FILTER Summary of this function goes here
%   Detailed explanation goes here
% input_args is the 1 dimensional signal
% ouput_args is the denoised 1 dimensional signal
H = form_Hankel( input_args );
[u s v]=svd(H, 'econ');

s1=zeros(size(s));
for i=1:size(s,1)
    s1(i,i)=(s(i,i)-r+abs(s(i,i)-r))/2;
end
H_denoised=u*s1*v';

sig_denoised=extract_1Dsignal_from_Hankel( H_denoised );% average 

end

