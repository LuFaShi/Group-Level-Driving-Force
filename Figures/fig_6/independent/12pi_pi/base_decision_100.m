% there are two subgroups with competing interests, termed as subgroup_1
% and subgroup_2, respectively. 

clear variables
clc

% the number of individuals
N=100;
% the proportions of informed individuals
p_1_sets=[0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08];
p_2=0.07;
% the weighting term for individuals
omega_1=0.3;
omega_2=0.3;
omega_3=0.3;
% the uncertainty range for the two informed subgroups
ran_1=pi/2;
ran_2=pi;
ran_3=2*pi; % uninformed individuals
% the preferred direction of informed individuals
g_1=[2 3];
g_1=g_1./sqrt(sum(g_1.^2));
% the angle between the preferred direction of subgroup_2 and the one of
% subgroup_1
theta_sets=3.0;
% the length of time step
delta_t=0.1;
% the time of evolution
T=500;
% the number of time step
Nt=T/delta_t;
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

% the frequency that the group selects g_1 or g_2 direction, respectively  
N_dir=zeros(length(p_1_sets),2);
% the mean effective magnitude of the group-level driving force
mag_force_sets=zeros(length(p_1_sets),2);

% the number of taking average
N_ave=1000;

% the number of computing worker
parpool(50);

for p_1_i=1:length(p_1_sets)
    p_1=p_1_sets(p_1_i);

    % the group direction
    h_sets=zeros(N_ave,length(theta_sets));
    diff_sets=zeros(N_ave,length(theta_sets));
    diff_sets_=zeros(N_ave,length(theta_sets));
    % the fragmentation probability of group
    fra_pro=zeros(N_ave,length(theta_sets));
    % the mean magnitude of the group-level driving force
    mag_force=zeros(N_ave,length(theta_sets));
    % the mean effective magnitude
    pro_1=zeros(N_ave,length(theta_sets));
    pro_2=zeros(N_ave,length(theta_sets));

    for theta_i=1:length(theta_sets)
        theta=theta_sets(theta_i);
        % the preferred direction of subgroup_2
        g_2_y=g_1(2)*cos(theta)+g_1(1)*sin(theta);
        g_2_x=g_1(1)*cos(theta)-g_1(2)*sin(theta);
        g_2=[g_2_x g_2_y];
        % the center of g_1 and g_2
        g_half_y=g_1(2)*cos(theta/2)+g_1(1)*sin(theta/2);
        g_half_x=g_1(1)*cos(theta/2)-g_1(2)*sin(theta/2);
        g_half=[g_half_x g_half_y];
     
        % the group direction
        h=zeros(N_ave,2);
        
        % counting the number that groups split
        fra=zeros(N_ave,1);
        % the mean magnitude of the group-level driving force in the ensemble
        mag_force_ensemble=zeros(N_ave,1);
      
        g_1_1=g_1(1);g_1_2=g_1(2);
        g_2_1=g_2(1);g_2_2=g_2(2);
        parfor ave_i=1:N_ave
            rng(2*ave_i-1,'twister')
            % the initial condition: random positions and orientations
            % the size of lattice
            s_lat=3;
            c_x=s_lat/2*(2*rand(N,1)-1);c_y=s_lat/2*(2*rand(N,1)-1);
            % the unit vector
            ran_angle=2*pi*rand(N,1);
            v_x=cos(ran_angle);v_y=sin(ran_angle);
            c=[c_x c_y]; v=[v_x v_y];
    
            % the centroid of the group
            cent=zeros(Nt,2);
    
            % the mean magnitude of the group-level driving force
            mag_force_=zeros(Nt,2);
            
            
            for t=1:Nt
               %% update the velocity of all individuals
                new_v=zeros(N,2); % all individuals require update simultaneously!
                sum_cos=0;
                sum_sin=0;
                for i=1:N
                    % define two summation vectors to calculate the velocity according to the
                    % repulsion and attraction rules, respectively
                    v_r=[0 0];
                    v_a=v(i,:)./sqrt(sum(v(i,:).^2));
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
                    % direction and their social interactions with weighting term
                    if i>=1 && i<=round(N*p_1)
                       theta_rot=(2*rand(1,1)-1)*ran_1/2; 
                       g_1_y_theta=g_1_2*cos(theta_rot)+g_1_1*sin(theta_rot);
                       g_1_x_theta=g_1_1*cos(theta_rot)-g_1_2*sin(theta_rot);
                       g_1_=[g_1_x_theta,g_1_y_theta];
                       
                       v_hat=(v_hat+omega_1*g_1_)./sqrt(sum((v_hat+omega_1*g_1_).^2));
                       sum_cos=sum_cos+g_1_(1);
                       sum_sin=sum_sin+g_1_(2);
                       
                    elseif i>=round(N*p_1)+1 && i<=round(N*p_1)+round(N*p_2)
                       theta_rot=(2*rand(1,1)-1)*ran_2/2; 
                       g_2_y_theta=g_2_2*cos(theta_rot)+g_2_1*sin(theta_rot);
                       g_2_x_theta=g_2_1*cos(theta_rot)-g_2_2*sin(theta_rot);
                       g_2_=[g_2_x_theta,g_2_y_theta];
                        
                        
                       v_hat=(v_hat+omega_2*g_2_)./sqrt(sum((v_hat+omega_2*g_2_).^2));
                       sum_cos=sum_cos+g_2_(1);
                       sum_sin=sum_sin+g_2_(2);
  
                    else
                       theta_rot=(2*rand(1,1)-1)*ran_3/2; 
                       g_half_y_theta=g_half_y*cos(theta_rot)+g_half_x*sin(theta_rot);
                       g_half_x_theta=g_half_x*cos(theta_rot)-g_half_y*sin(theta_rot);
                       g_half_=[g_half_x_theta,g_half_y_theta];
                         
                       v_hat=(v_hat+omega_3*g_half_)./sqrt(sum((v_hat+omega_3*g_half_).^2)); 
                       sum_cos=sum_cos+g_half_(1);
                       sum_sin=sum_sin+g_half_(2);
                       
                    end
    
                    % rotating the velocity vector by a random angle because the
                    % individual is subject to random influences
                    % the random angle is drawn from the Gaussian distribution, centred
                    % on 0, with standard deviation sigma radians--sqrt(sigma)*randn
    
                    % In three dimensional space, we firstly need to determine the 
                    % plane of rotating, and we assume that the probabilities of
                    % rotating are same in any plane. The plane is represented by the velocity
                    % vector itself and a vector located on x-y plane, referred as vec_xy
    
                    % the angle of rotating
                    theta_rot_2=sqrt(sigma)*randn(1,1);
                    v_hat_y=v_hat(2)*cos(theta_rot_2)+v_hat(1)*sin(theta_rot_2);
                    v_hat_x=v_hat(1)*cos(theta_rot_2)-v_hat(2)*sin(theta_rot_2);
                    v_hat=[v_hat_x v_hat_y];
    
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
                % the group-level driving force
                mag_force_(t,:)=[sum_cos sum_sin]./N;
                   
            end
            % examine whether the group split
            for kk=2:N
               if sqrt(sum((c(kk,:)-c(1,:)).^2))>10*rho
                    fra(ave_i)=fra(ave_i)+1;
                    break
               end
            end
            h(ave_i,:)=(cent(Nt,:)-cent(Nt-100,:))./sqrt(sum((cent(Nt,:)-cent(Nt-100,:)).^2));

            % the projection of the group-level driving force on the two
            % preferred directions, i.e., the mean effective magnitude
            pro_1(ave_i)=mean(sum(mag_force_.*g_1,2));
            pro_2(ave_i)=mean(sum(mag_force_.*g_2,2));
            
        end
        
        h_sets(:,theta_i)=atan2(h(:,2),h(:,1));
        diff_sets(:,theta_i)=acos(sum(h.*g_1,2)./(sqrt(sum(h.^2,2)).*sqrt(sum(g_1.^2))));
        diff_sets_(:,theta_i)=acos(sum(h.*g_2,2)./(sqrt(sum(h.^2,2)).*sqrt(sum(g_2.^2))));
        fra_pro(:,theta_i)=fra;
    end

    % counting the number that the group selects a specific direction
    s1=0;s2=0;
    F_1=0;F_2=0;
    for ii=1:N_ave
        if fra_pro(ii)~=1
            if diff_sets(ii,length(theta_sets))<diff_sets_(ii,length(theta_sets))
                s1=s1+1;
                F_1=F_1+pro_1(ii);
            else
                s2=s2+1;
                F_2=F_2+pro_2(ii);
            end
        end
    end

    N_dir(p_1_i,:)=[s1 s2];
    mag_force_sets(p_1_i,:)=[F_1 F_2];

end
save('please enter the file path','-regexp','N_dir','mag_force_sets','fra_pro')
