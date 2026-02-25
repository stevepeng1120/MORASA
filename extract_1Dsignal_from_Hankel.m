function [ output_args ] = extract_1Dsignal_from_Hankel( input_args )
% EXTRACT_1DSIGNAL_FROM_HANKEL Summary of this function goes here
% Detailed explanation goes here
% Author: Xi Peng
H=input_args;
[Nx,Ny]=size(H);
signal_length=Nx+Ny-1;
signal_matrix=zeros(Nx,signal_length);% row vector
for i=1:Nx
    signal_matrix(i,i:i+Ny-1) = H(i,:);
end

if mod(signal_length,2)==0
    ave_weight=[1:1:signal_length/2, signal_length/2:-1:1];
end
if mod(signal_length,2)==1
    ave_weight=[1:1:(signal_length+1)/2, (signal_length-1)/2:-1:1];
end
output_args=sum(signal_matrix,1)./ave_weight;
output_args=output_args(:);

end

