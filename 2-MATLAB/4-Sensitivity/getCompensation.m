function [ V ] = getCompensation( P, i, j, ctype )
%GETCOMPENSATION returns the compensation pattern matrix
%   used with Pstat_deriv
%
%INPUT
% P - matrix of transition probabilities
% i,j - the P element being perturbed primarily
% ctype - compensation pattern to maintain column stochasticity of P
%       1 - Specific: dp_ij = 1; dp_jj = -1
%       2 - Uniform: dp_ij = 1; dp_mj = - 1 / (s-1)
%       3 - Proportional
%
% Output
% V - a matrix with zeros everywhere except at v_ij and one or some elements in column j

s = size(P,1);
V = zeros(s);

if ctype == 1                   %Specific
    V(i,j) = 1;
    V(j,j) = -1;
elseif ctype == 2               %Uniform
    for k=1:s
        V(k,j) = -1 / (s-1);
    end
    V(i,j) = 1;
else                            %Proportional
    for k=1:s
        V(k,j) = - P(k,j) / (1 - P(i,j));
    end
    V(i,j) = 1;
end

end

