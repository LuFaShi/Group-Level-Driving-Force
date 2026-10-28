clear variables
clc

% the number of individuals
N_sets=10;
% the weighting term for individuals
omega=0.5;
% the size of the uncertainty range within the two-component group, the
% control parameter is denoted by p, the proportion of informed individuals
p_sets=0.1:0.1:1.0;
% the preferred direction of individuals
g=[2 0];
g=g./sqrt(sum(g.^2));
% the length of time step
delta_t=0.1;
% the time of evolution
T=500;
% the number of time step
Nt=T/delta_t;
% the number of transient time
N_transient=500;
% the minimum distance between two individuals (range of repulsion)
alpha=1;
% the local interaction range
rho=6;
% the speed 
s=1;
% the turning rate
tur_rat=2.0;
% a threshold for updating velocity
thre=delta_t*tur_rat;
% the variance of the Gaussian distribution
sigma=0.01;

% the number of intervals for the statistical distribution of deviations 
% from the preferred direction
N_int=20;
% the number of taking ensemble average
N_ensemble=1000;

% We adopt two measures to quantify the performance of information
% transfer, one is the accuracy used in the original paper, another is the
% order parameter
acc_sets=zeros(length(N_sets),length(p_sets));
OP_sets=zeros(length(N_sets),length(p_sets));
% the proportion that groups split (fragment)
pro_spl=zeros(length(N_sets),length(p_sets));
% the probability density of deviations
pro_den=zeros(N_int,length(N_sets),length(p_sets));
% the mean value
mea_sets=zeros(length(N_sets),length(p_sets));
% the variance
vari_sets=zeros(length(N_sets),length(p_sets));
% the three-order original moment
three_ori_moment_sets=zeros(length(N_sets),length(p_sets));
% the three-order center moment
three_center_moment_sets=zeros(length(N_sets),length(p_sets));
% the skewness
skewness_sets=zeros(length(N_sets),length(p_sets));

% the group-level driving force
mag_sets=zeros(length(N_sets),length(p_sets));
mag_var_sets=zeros(length(N_sets),length(p_sets));
dir_sets=zeros(length(N_sets),length(p_sets));
dir_var_sets=zeros(length(N_sets),length(p_sets));
% the center-of-mass velocity
mag_cent_sets=zeros(length(N_sets),length(p_sets));
mag_var_cent_sets=zeros(length(N_sets),length(p_sets));
dir_cent_sets=zeros(length(N_sets),length(p_sets));
dir_var_cent_sets=zeros(length(N_sets),length(p_sets));

% the number of computing worker
parpool(50);

for N_i=1:length(N_sets)
    N=N_sets(N_i);
    for p_i=1:length(p_sets)       
        p=p_sets(p_i);

        %% physical quantity in the ensemble
        % the group direction 
        h=zeros(N_ensemble,2);
        % the order parameter
        OP_ens=zeros(N_ensemble,1);
        % the probability distribution of the deviation
        P_ens=zeros(N_ensemble,N_int);
        % the mean value
        m_ens=zeros(N_ensemble,1);
        % the variance
        var_ens=zeros(N_ensemble,1);
        % the three-order moment
        three_ori_moment_ens=zeros(N_ensemble,1);
        three_center_moment_ens=zeros(N_ensemble,1);
        skewness_ens=zeros(N_ensemble,1);
        % counting the number that groups split
        Index_split=zeros(N_ensemble,1);
        % the group-level driving force
        mag_ens=zeros(N_ensemble,1);
        mag_var_ens=zeros(N_ensemble,1);
        dir_ens=zeros(N_ensemble,1);
        dir_var_ens=zeros(N_ensemble,1);
        % the center-of-mass velocity
        mag_cent_ens=zeros(N_ensemble,1);
        mag_var_cent_ens=zeros(N_ensemble,1);
        dir_cent_ens=zeros(N_ensemble,1);
        dir_var_cent_ens=zeros(N_ensemble,1);

        %% ************ end **********

        % g is a broadcast variable
        g_1=g(1);g_2=g(2);
        parfor ensemble_i=1:N_ensemble
            rng(2*ensemble_i-1,'twister')
            % the initial condition: random positions and orientations
            s_lat=3;
            c_x=s_lat/2*(2*rand(N,1)-1);c_y=s_lat/2*(2*rand(N,1)-1);
            % the unit vector
            ran_angle=2*pi*rand(N,1);
            v_x=cos(ran_angle);v_y=sin(ran_angle);
            c=[c_x c_y]; v=[v_x v_y];
            % the uncertainty range
            ran=2*pi*ones(N,1);
            ran(1:round(N*p))=0;

            % the Order Parameter
            OP=zeros(Nt,1);
            % the individuals' deviations from the preferred direction 
            dev=zeros(Nt*N,1);
            % the center-of-mass the group
            cent=zeros(Nt,2);
            % the group-level driving force
            mag=zeros(Nt,1);
            dir=zeros(Nt,1);
            % the center-of-mass velocity
            mag_cent=zeros(Nt,1);
            dir_cent=zeros(Nt,1);
            for t=1:Nt

               new_v=zeros(N,2); % all individuals require update simultaneously!
               force_x=zeros(1,1);
               force_y=zeros(1,1);

               %% update the velocity of all individuals 
                for i=1:N
                    % define two summation vectors to calculate the velocity according to the
                    % repulsion and attraction rules, respectively
                    v_r=[0 0];
                    v_a=v(i,:);
                    for j=1:N
                        % the distance between the focal individual and the other
                        dis=sqrt(sum((c(i,:)-c(j,:)).^2));
                        if dis>0 && dis<=alpha
                           v_r=v_r-(c(j,:)-c(i,:))./sqrt(sum((c(j,:)-c(i,:)).^2));
                        elseif dis>alpha && dis<=rho
                           v_a=v_a+(c(j,:)-c(i,:))./sqrt(sum((c(j,:)-c(i,:)).^2))+...
                               v(j,:)./sqrt(sum(v(j,:).^2));
                        end            
                    end
                    % determine the interaction rule applied by the focal individual 
                    if sum(v_r.^2)>0
                        v_hat=v_r./sqrt(sum(v_r.^2));
                    else
                        v_hat=v_a./sqrt(sum(v_a.^2));
                    end

                    % individuals balance the influence of their preferred
                    % direction and their social interactions with a
                    % weighting factor, the dependent scheme

                    theta_=acos(sum(v_hat.*g)/(sqrt(sum(v_hat.^2))*sqrt(sum(g.^2))));
                    if theta_<ran(i)/2
                        g_=v_hat;
                        dev((t-1)*N+i)=theta_;
                    else
                        g_=rotating([g 0],[v_hat 0],ran(i)/2);
                        g_=g_(1:2);
                        dev((t-1)*N+i)=ran(i)/2;
                    end
                    force_x=force_x+g_(1);
                    force_y=force_y+g_(2);

                    v_hat=(v_hat+omega*g_)./sqrt(sum((v_hat+omega*g_).^2)); 

                    % rotating the velocity vector by a random angle because the
                    % individual is subject to random influences
                    % the random angle is drawn from the Gaussian distribution, centred
                    % on 0, with standard deviation sigma radians--sqrt(sigma)*randn

                    % the angle of rotating
                    theta_rot_2=sqrt(sigma)*randn(1,1);
                    % using the defined function, rotating
                    v_hat_y=v_hat(2)*cos(theta_rot_2)+v_hat(1)*sin(theta_rot_2);
                    v_hat_x=v_hat(1)*cos(theta_rot_2)-v_hat(2)*sin(theta_rot_2);
                    v_hat=[v_hat_x,v_hat_y];

                    % whether the difference between the updated velocity and the previous
                    % one exceeds the thre, if not, the individual achieves
                    % alignment with the updated velocity, otherwise the individual
                    % turn thre towards it
                    if acos(sum(v_hat.*v(i,:))/(sqrt(sum(v_hat.^2))*sqrt(sum(v(i,:).^2))))<=thre
                        new_v(i,:)=v_hat;
                    else
                        new_v_=rotating([v(i,:) 0],[v_hat 0],thre);
                        new_v(i,:)=new_v_(1:2);
                    end

                end

               %% update the position of all individuals
                 c=c+new_v*s*delta_t;
                 v=new_v;

                 cent(t,:)=sum(c)/N;
                 % the order parameter
                 OP(t)=sqrt(sum(mean(v(:,1))^2+mean(v(:,2))^2));
                 % the group-level driving force
                 mag(t)=force_x/N;
                 dir(t)=force_y/N;
                 % the center-of-mass velocity
                 mag_cent(t)=mean(v(:,1));
                 dir_cent(t)=mean(v(:,2));
                 if t>N_transient
                    % examine whether the group split
                    dis_pairs=zeros((N-1)*N/2,1);
                     for ii=1:N-1
                         for jj=ii+1:N
                             dis_pairs((ii-1)*N-((ii-1)*ii)/2+jj-ii)=sqrt(sum((c(ii,:)-c(jj,:)).^2));
                         end
                     end
                    for kk=1:length(dis_pairs)
                       if dis_pairs(kk)>(N*alpha)*2
                           Index_split(ensemble_i)=1;
                           break
                       end
                    end
                 end

            end

            % the group direction
            h(ensemble_i,:)=(cent(Nt,:)-cent(Nt-100,:))./sqrt(sum((cent(Nt,:)-cent(Nt-100,:)).^2));
            OP_ens(ensemble_i)=mean(OP(N_transient:end));
            % the statistical distribution of deviations
            [N_set,P_set,x_sca]=sta_dis(abs(dev),N_int,pi);
            P_ens(ensemble_i,:)=P_set;
            % the mean deviation
            m_ens(ensemble_i)=mean(abs(dev));
            % the variance
            var_ens(ensemble_i)=mean((abs(dev)).^2)-mean(abs(dev))^2;
            % the three-order moment
            three_ori_moment_ens(ensemble_i)=mean(abs(dev).^3);
            three_center_moment_ens(ensemble_i)=mean((abs(dev)-mean(abs(dev))).^3);
            skewness_ens(ensemble_i)=three_center_moment_ens(ensemble_i)/(sqrt(var_ens(ensemble_i))^3);
            % the group-level driving force
            mag_ens(ensemble_i)=mean(mag(N_transient:end));
            mag_var_ens(ensemble_i)=mean(mag(N_transient:end).^2)-mean(mag(N_transient:end))^2;
            dir_ens(ensemble_i)=mean(dir(N_transient:end));
            dir_var_ens(ensemble_i)=mean(dir(N_transient:end).^2)-mean(dir(N_transient:end))^2;
            % the center-of-mass velocity
            mag_cent_ens(ensemble_i)=mean(mag_cent(N_transient:end));
            mag_var_cent_ens(ensemble_i)=mean(mag_cent(N_transient:end).^2)-mean(mag_cent(N_transient:end))^2;
            dir_cent_ens(ensemble_i)=mean(dir_cent(N_transient:end));
            dir_var_cent_ens(ensemble_i)=mean(dir_cent(N_transient:end).^2)-mean(dir_cent(N_transient:end))^2;

        end
        Index_nosplit=find(Index_split==0);
        S=N_ensemble-length(Index_nosplit);
        % the accuracy
        % the difference between the group direction and the preferred direction
        diff_g=acos(sum(h(Index_nosplit,:).*g,2)./(sqrt(sum(h(Index_nosplit,:).^2,2)).*sqrt(sum(g.^2))));
        acc_sets(N_i,p_i)=1-mean(diff_g)/(pi/2);
        % the order parameter
        OP_sets(N_i,p_i)=mean(OP_ens(Index_nosplit));
        % the proportion
        pro_spl(N_i,p_i)=S/N_ensemble;
        % the probability density
        pro_den(:,N_i,p_i)=mean(P_ens(Index_nosplit,:))';
        % the mean value
        mea_sets(N_i,p_i)=mean(m_ens(Index_nosplit));
        % the variance
        vari_sets(N_i,p_i)=mean(var_ens(Index_nosplit));
        % the three-order moment
        three_center_moment_sets(N_i,p_i)=mean(three_center_moment_ens(Index_nosplit));
        three_ori_moment_sets(N_i,p_i)=mean(three_ori_moment_ens(Index_nosplit));
        skewness_sets(N_i,p_i)=mean(skewness_ens(Index_nosplit));
        % the group-level driving force
        mag_sets(N_i,p_i)=mean(mag_ens(Index_nosplit));
        mag_var_sets(N_i,p_i)=mean(mag_var_ens(Index_nosplit));
        dir_sets(N_i,p_i)=mean(dir_ens(Index_nosplit));
        dir_var_sets(N_i,p_i)=mean(dir_var_ens(Index_nosplit));
        % the center-of-mass velocity
        mag_cent_sets(N_i,p_i)=mean(mag_cent_ens(Index_nosplit));
        mag_var_cent_sets(N_i,p_i)=mean(mag_var_cent_ens(Index_nosplit));
        dir_cent_sets(N_i,p_i)=mean(dir_cent_ens(Index_nosplit));
        dir_var_cent_sets(N_i,p_i)=mean(dir_var_cent_ens(Index_nosplit));

    end

end
save('/share/home/shilufa/couzin/phenomenon_knowledge/dependent/omega/0_5/two_component/acc.mat','-regexp',...
    'acc_sets','OP_sets','pro_spl','pro_den','mea_sets','vari_sets',...
    'mag_sets','mag_var_sets','dir_sets','dir_var_sets',...
    'mag_cent_sets','mag_var_cent_sets','dir_cent_sets','dir_var_cent_sets',...
    'three_center_moment_sets','three_ori_moment_sets','skewness_sets');

