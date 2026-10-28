% this function is used to calculate the statistical distribution of
% deviations from the preferred direction

% the input: 
% a random variable (ran_var)
% the number of intervals (N_int)
% the maximum variable range (max_ran) 

% the output:
% the number of the random variable appeared in each interval,called N_int; 
% the corresponding probability, called P_int
% the scale value of x axis, called x_sca
function [N_set,P_set,x_sca]=sta_dis(ran_var,N_int,max_ran)

        % the scale of intervals
        sca_int=linspace(0,max_ran,N_int+1);
        % the length of intervals
        len_int=max_ran/N_int; 
        % define a vector to store the number of the random variable 
        % appeared in each interval
        N_set=zeros(N_int,1);
        for i=1:length(ran_var)
            for j=1:N_int
                if ran_var(i,1)>=sca_int(j) && ran_var(i,1)<sca_int(j+1)
                    N_set(j)=N_set(j)+1;
                    break
                end
            end    
        end
        % the probability
        P_set=N_set./sum(N_set);
        P_set=P_set./len_int;
        % the x_scale
        x_sca=sca_int(2:end)-len_int/2;
end