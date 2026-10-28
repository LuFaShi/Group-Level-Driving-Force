clear variables
clc

% the number of individuals
N=50;
% the weighting factor for informed individuals
omega=0.4;
% the proportion of informed individuals
p_sets=(1:3:N)./N;
% the preferred direction of informed individuals
g=[-2 -3];
g=g./sqrt(sum(g.^2));
% the length of time step
delta_t=0.1;
% the time of evolution
T=500;
% the number of time step
Nt=T/delta_t;
% the minimum distance between two individuals (range of repulsion)
alpha=1;
% the local interaction range
rho_sets=5:0.5:8;
% the speed 
s=1.0;
% the turning rate
tur_rat=2;
% a threshold for updating velocity
thre=delta_t*tur_rat;
% the variance of the Gaussian distribution
sigma=0.01;
% the number of taking ensemble average
N_ensemble=500;

% fluctuations of directions of center-of-mass velocity and Langevin-like force
flu_vel=zeros(length(rho_sets),length(p_sets));
flu_force=zeros(length(rho_sets),length(p_sets));
% the accuracy of collective decision-making
acc=zeros(length(rho_sets),length(p_sets));
% the proportion of group split (fragmentation)
pro_spl=zeros(length(rho_sets),length(p_sets));
% the average neighbour number of individuals
ave_N_nei=zeros(length(rho_sets),length(p_sets));
% the number of computing worker
parpool(20)
g_2=g(2);g_1=g(1); % In parallel computing, g is a broadcast variable
for rho_i=1:length(rho_sets)
    rho=rho_sets(rho_i);
    for p_i=1:length(p_sets)
        p=p_sets(p_i);

        % physical quantity in the ensemble
        flu_vel_ens=zeros(N_ensemble,1);
        flu_force_ens=zeros(N_ensemble,1);
        h=zeros(N_ensemble,2);
        ave_N_nei_ens=zeros(N_ensemble,1);
        S=0; % count number of group splitting

        parfor ensemble_i=1:N_ensemble
            rng(2*ensemble_i-1,'twister')

            % the initial condition: random positions and directions
            % the size of lattice
            s_lat=3;
            c_x=s_lat/2*(2*rand(N,1)-1);c_y=s_lat/2*(2*rand(N,1)-1);
            % the unit vector
            ran_angle=2*pi*rand(N,1);
            v_x=cos(ran_angle);v_y=sin(ran_angle);
            c=[c_x c_y]; v=[v_x v_y];
            
            % the knowledge level of individuals is quantified by the szie of an
            % uncertain angular range centered on the accurate direction, representing
            % the desired traveling route
            % the larger the size of the uncertain range, the lower the knowledge level
            
            % here, the uncertain range of all individuals is determined by 
            % the proportion of informed individuals, p 

            % the range
            ran=2*pi*ones(N,1);
            ran(1:round(N*p))=0;
            % the center-of-mass of the group
            cent=zeros(Nt,2);
            % the deviation of actual preferences of all individuals from the desired one 
            dev=zeros(Nt*N,1);
            
            % the social interactions
            soc_int=zeros(Nt*N,2);
            % the Langevin-like force
            Lan_force=zeros(Nt,2);
            % the center-of-mass velocity
            cen_vel=zeros(Nt,2);
            % count the number of interacting occurred in each time step 
            N_int=zeros(Nt,1);
            for t=1:Nt
            
                %% update the velocity of all individuals
                new_v=zeros(N,2); % all individuals require update simultaneously!
            
                Lan_force_x=0;
                Lan_force_y=0;
                
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
                           N_int(t)=N_int(t)+1;
                        end            
                    end
                    % determine the social interaction acting on the focal individual 
                    if sum(v_r.^2)>0
                        v_hat=v_r./sqrt(sum(v_r.^2));
                    else
                        v_hat=v_a./sqrt(sum(v_a.^2));
                    end
                    % recording the social interactions
                    soc_int((t-1)*N+i,:)=v_hat;
            
                    % individuals balance the influence of their preferred
                    % direction and their social interactions with a weighting factor
                    % the independent scheme
                    
                    theta_rot_1=(2*rand(1,1)-1)*(ran(i)/2);
                    dev((t-1)*N+i)=theta_rot_1;
             
                    g_y=g_2*cos(theta_rot_1)+g_1*sin(theta_rot_1);
                    g_x=g_1*cos(theta_rot_1)-g_2*sin(theta_rot_1);
                    g_=[g_x g_y];
                    Lan_force_x=Lan_force_x+g_x;
                    Lan_force_y=Lan_force_y+g_y;
            
                    v_hat=(v_hat+omega*g_)./sqrt(sum((v_hat+omega*g_).^2)); 
                    
                    % rotating the velocity vector by a random angle because the
                    % individual is subject to random errors
                    % the random angle is drawn from the Gaussian distribution, centred
                    % on 0, with standard deviation sigma radians--sqrt(sigma)*randn
            
                    % the angle of rotating
                    theta_rot_2=sqrt(sigma)*randn(1,1);
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
                % the center-of-mass position 
                cent(t,:)=sum(c)/N;
                % the center-of-mass velocity and the Langevin-like force 
                cen_vel(t,:)=[mean(v(:,1)),mean(v(:,2))];
                Lan_force(t,:)=[Lan_force_x/N,Lan_force_y/N];            
             
            end
            % the directions of the center-of-mass velocity and the
            % Langevin-like force
            d_cen_vel=atan2(cen_vel(:,2),cen_vel(:,1));
            d_Lan_force=atan2(Lan_force(:,2),Lan_force(:,1));
            d_cen_vel=mod(d_cen_vel-atan2(g_2,g_1)+pi,2*pi);
            d_Lan_force=mod(d_Lan_force-atan2(g_2,g_1)+pi,2*pi);
            flu_force_ens(ensemble_i)=var(d_Lan_force);
            flu_vel_ens(ensemble_i)=var(d_cen_vel);
            % the direction of group movement
            h(ensemble_i,:)=(cent(Nt,:)-cent(Nt-100,:))./sqrt(sum((cent(Nt,:)-cent(Nt-100,:)).^2));
            % examine whether the group split
            for kk=2:N
               if sqrt(sum((c(kk,:)-c(1,:)).^2))>10*rho
                   S=S+1;
                   break
               end
            end
            % the average number of interacting neighbours of each dividual
            % at each time step after discarding transient
            ave_N_nei_ens(ensemble_i)=mean(N_int(100:end)./N);

        end

        flu_vel(rho_i,p_i)=mean(flu_vel_ens);
        flu_force(rho_i,p_i)=mean(flu_force_ens);
        pro_spl(rho_i,p_i)=S/N_ensemble;
        % the angular difference between the group direction and the
        % preferred direction
        ang_diff=acos(sum(h.*g,2)./(sqrt(sum(h.^2,2))*sqrt(sum(g.^2))));
        acc(rho_i,p_i)=1-mean(ang_diff)/(pi/2);
        ave_N_nei(rho_i,p_i)=mean(ave_N_nei_ens);
    end
end
save('/share/home/shilufa/couzin_/dependence_rho/independent/acc.mat','-regexp',...
    'flu_vel','flu_force','pro_spl','acc','ave_N_nei')
