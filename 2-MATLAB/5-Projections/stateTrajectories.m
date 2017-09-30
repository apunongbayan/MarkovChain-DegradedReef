function [ trajectories, tvec ] = stateTrajectories( Psamples, x0, time_fin )
%MC_TRAJECTORIES - Projections of state covers given initial condition
%Input
%   Psamples - sampled transition probabilities (row - sample, col - pij)
%   x0 - vector N by 1 of initial state distribution
%   time_fin - final time point; in months
%Output
%   trajectories - 3d matrix; x - state, y - time, z - replicate
%   tvec - time points

a = size(x0);                %check that it is a col vector
if a(2) > 1
    x0 = x0';                %transpose
end

[nr, nc] = size(Psamples); 

trajectories = [];          % 3d matrix

for i = 1:nr
    trajectory = [x0 ];
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

