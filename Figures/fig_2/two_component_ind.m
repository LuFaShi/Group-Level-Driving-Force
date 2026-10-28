clear variables
clc

% the number of individuals
N_sets=[10 30 50 100];
% the weighting term for individuals
omega=0.4;
% the size of the uncertainty range within the two-component group, the
% control parameter is denoted by p, the proportion of informed individuals
% Since the proportion in the two-component group is different for
% different group sizes, we set a maximum matrix to contain all elements
% under arbitrary group size for simplicity
p_sets=(1:5:100)./100;
% the preferred direction of individuals
g=[2 3];
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
tur_rat=2;
% a threshold for updating velocity
thre=delta_t*tur_rat;
% the variance of the Gaussian distribution
sigma=0.01;

% the number of intervals for the statistical distribution of deviations 
% from the preferred direction
N_int=20;
% the number of taking ensemble average
N_ensemble=2000;

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

% the number of computing worker
parpool(50);

for N_i=1:length(N_sets)
    N=N_sets(N_i);
    if N==10
        p_sets_prime=(1:1:10)./10;
    elseif N==30
        p_sets_prime=(1:2:30)./30;
    elseif N==50
        p_sets_prime=(1:3:50)./50;
    else
        p_sets_prime=(1:5:100)./100;
    end
    for p_i=1:length(p_sets_prime)       
        p=p_sets_prime(p_i);

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
        % counting the number that groups split
        Index_split=zeros(N_ensemble,1);
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
            for t=1:Nt

               new_v=zeros(N,2); % all individuals require update simultaneously!

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
                    % weighting factor, the independent scheme

                    theta_rot_1=(2*rand(1,1)-1)*(ran(i)/2);
                    dev((t-1)*N+i)=theta_rot_1;
                    g_y=g_2*cos(theta_rot_1)+g_1*sin(theta_rot_1);
                    g_x=g_1*cos(theta_rot_1)-g_2*sin(theta_rot_1);
                    g_=[g_x g_y];
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
    end

end
save('please enter the file path','-regexp',...
    'acc_sets','OP_sets','pro_spl','pro_den','mea_sets','vari_sets')
