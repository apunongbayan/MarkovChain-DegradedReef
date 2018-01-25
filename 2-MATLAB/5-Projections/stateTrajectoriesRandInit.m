function [ trajectories, tvec ] = stateTrajectoriesRandInit( Psamples, Xsamples, time_fin )
%STATETRAJECTORIESRANDINIT Projections of state covers, random initial
%distributions
%   Input
%       Psamples - sampled transition probabilities (row - sample, col - pij)
%       Xsamples - sampled initial state distributions (row - sample, col -
%       initial cover of state j)         
%       time_fin - final time point; in months
%   Output
%       trajectories - 3d matrix; x - state, y - time, z - replicate
%       tvec - time points; in months

[nr, nc] = size(Psamples); 
trajectories = [];          % 3d matrix

for i = 1:nr
    trajectory = [Xsamples(i,:)' ];
    tvec = [0 ];
    P = reshape(Psamples(i,:),8,8);
    for k = 2:time_fin+1
        state = P * trajectory(:,k-1);
        trajectory = [trajectory state];
        tvec(end+1) = k-1;
    end
    trajectories(:,:,i) = trajectory;
end



end

