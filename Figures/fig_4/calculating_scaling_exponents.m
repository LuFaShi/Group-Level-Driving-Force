% This code is used to calculate the scaling exponents of the y-direction
% variance of both the group-level driving force and the center-of-mass
% velocity; "b_sets" is the control parameter sets of the homogeneous group, 
% it can be replaced by "k_sets" corresponding to the powerlaw group or other groups

% calculating the slope
slope_sets=zeros(length(b_sets),1);
% the fitting goodness
goodness=zeros(length(b_sets),1);

for ii=1:length(b_sets)
    [a,b]=linear_reg(log(N_sets),log(dir_var_cent_sets(:,ii)));
    slope_sets(ii)=b;
    % the predicted value
    pre_value=log(N_sets).*b+a;
    % the difference between the fitting value and the actual one
    diff_1=sum((pre_value'-log(dir_var_cent_sets(:,ii))).^2);
    % the difference between the actual value and its mean
    diff_2=sum((log(dir_var_cent_sets(:,ii))-mean(log(dir_var_cent_sets(:,ii)))).^2);
    % the fitting goodness
    goodness(ii)=1-(diff_1/diff_2);
end
plot(mea_sets(1,:),abs(slope_sets),'-o')
hold on
plot(mea_sets(1,:),goodness,'-o')