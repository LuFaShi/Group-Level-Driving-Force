clear variables
clc

% the number of individuals
N=50;
% the weighting term for individuals
omega=0.4;
% the preferred direction of individuals
g=[2 3];
g=g./sqrt(sum(g.^2));
% the length of time step
delta_t=0.1;
% the time of evolution
T=5000;
% the number of time step
Nt=T/delta_t;
% the minimum distance between two individuals (range of repulsion)
alpha=1;
% the local interaction range
rho=6;
% the speed 
s=1.0;
% the turning rate
tur_rat=2;
% a threshold for updating velocity
thre=delta_t*tur_rat;
% the variance of the Gaussian distribution
sigma=0.01;

rng(1,'twister')
% the initial condition: random positions and directions
% the size of lattice
s_lat=3;
c_x=s_lat/2*(2*rand(N,1)-1);c_y=s_lat/2*(2*rand(N,1)-1);
% the unit vector
ran_angle=2*pi*rand(N,1);
v_x=cos(ran_angle);v_y=sin(ran_angle);
c=[c_x c_y]; v=[v_x v_y];
 
% the proportion of informed individuals 
p=0.1;
% the uncertainty range
ran=2*pi*ones(N,1);
ran(1:round(N*p))=0;
% the center-of-mass of the group
cent=zeros(Nt,2);
% the individuals' deviations from the preferred direction 
dev=zeros(Nt*N,1);

% the social interactions
soc_int=zeros(Nt*N,2);
% the group-level driving force
group_force=zeros(Nt,2);
% the center-of-mass velocity
cen_vel=zeros(Nt,2);

for t=1:Nt

    %% update the velocity of all individuals
    new_v=zeros(N,2); % all individuals require update simultaneously!

    force_x=0;
    force_y=0;
    
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
        % recording the social interactions
        soc_int((t-1)*N+i,:)=v_hat;

        % individuals balance the influence of their preferred
        % direction and their social interactions with a
        % weighting factor, the independent scheme

        theta_rot_1=(2*rand(1,1)-1)*(ran(i)/2);
        dev((t-1)*N+i)=theta_rot_1;
        g_y=g(2)*cos(theta_rot_1)+g(1)*sin(theta_rot_1);
        g_x=g(1)*cos(theta_rot_1)-g(2)*sin(theta_rot_1);
        force_x=force_x+g_x;
        force_y=force_y+g_y;
        g_=[g_x g_y];
        v_hat=(v_hat+omega*g_)./sqrt(sum((v_hat+omega*g_).^2));
        
        % rotating the velocity vector by a random angle because the
        % individual is subject to random errors
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
    % the center-of-mass position 
    cent(t,:)=sum(c)/N;
    % the center-of-mass velocity and the group-level-driving force 
    cen_vel(t,:)=[mean(v(:,1)),mean(v(:,2))];
    group_force(t,:)=[force_x/N,force_y/N];

 
end
% the directions of the center-of-mass velocity, the group-level-driving force
% and the social interactions
d_cen_vel=atan2(cen_vel(:,2),cen_vel(:,1));
d_group_force=atan2(group_force(:,2),group_force(:,1));
d_soc_int=atan2(soc_int(:,2),soc_int(:,1));

% taking modulus of 2*pi, for the sake of clarity of the plot
d_group_force=mod(d_group_force-atan2(g(2),g(1))+pi,2*pi);
d_cen_vel=mod(d_cen_vel-atan2(g(2),g(1))+pi,2*pi);
d_soc_int=mod(d_soc_int-atan2(g(2),g(1))+pi,2*pi);
% 
% using the function(sta_dis) to calculate the distribution
figure(1)
[~,P_set,x_sca]=sta_dis(d_cen_vel,80,2*pi);
plot(x_sca,P_set,'o')
hold on
[~,P_set,x_sca]=sta_dis(d_group_force,80,2*pi);
plot(x_sca,P_set,'o')

figure(2)
[~,P_set,x_sca]=sta_dis(d_soc_int,80,2*pi);
plot(x_sca,P_set,'o')

figure(3)
plot((1:T/delta_t)*delta_t,d_cen_vel)
hold on
plot((1:T/delta_t)*delta_t,d_group_force)

figure(4)
plot(cent(:,1),cent(:,2),'LineWidth',2)

