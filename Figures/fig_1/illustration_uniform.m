clear variables
clc

% the number of individuals
N=30;
% the weighting term for individuals
omega=0.4;
% the preferred direction of informed individuals
g=[2 3];
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
rho=6;
% the speed 
s=1;
% the turning rate
tur_rat=2;
% a threshold for updating velocity
thre=delta_t*tur_rat;
% the variance of the Gaussian distribution
sigma=0.01;

% rng(1,'twister')
% the initial condition: random positions and orientations
s_lat=3;
c_x=s_lat/2*(2*rand(N,1)-1);c_y=s_lat/2*(2*rand(N,1)-1);
% the unit vector
ran_angle=2*pi*rand(N,1);
v_x=cos(ran_angle);v_y=sin(ran_angle);
c=[c_x c_y]; v=[v_x v_y];

% the uncertainty range; a uniform distribution with range of 2*pi
ran=(2*pi).*rand(N,1);

% the center-of-mass
cent=zeros(Nt,2);

% the deviation of an individual's actual preferred direction from the
% accurate direction
dev=zeros(Nt*N,1);
for t=1:Nt

    new_v=zeros(N,2); % all individuals require update simultaneously!
%      quiver(c(:,1),c(:,2),v(:,1),v(:,2),0.3);
%      pause(0.1)
    
   %% update the velocity of all individuals 
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
        % direction and their social interactions with a
        % weighting factor, the independent scheme

        theta_rot_1=(2*rand(1,1)-1)*(ran(i)/2);        
        dev((t-1)*N+i)=theta_rot_1;
        g_y=g(2)*cos(theta_rot_1)+g(1)*sin(theta_rot_1);
        g_x=g(1)*cos(theta_rot_1)-g(2)*sin(theta_rot_1);
        g_=[g_x g_y];
        v_hat=(v_hat+omega*g_)./sqrt(sum((v_hat+omega*g_).^2)); 
        
        % rotating the velocity vector by a random angle because the
        % individual is subject to random influences
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
     
    cent(t,:)=sum(c)/N;
 

end
colormap("jet")
colorlist=jet;
% the extent of knowledge is divided according to colorlist
ran_ext=linspace(0,2*pi,size(colorlist,1));
for p_i=1:N
    % the difference between the actual extent and the certain discrete value
    diff_ran=abs(ran(p_i)-ran_ext);
    % the index corresponding to the specific color
    [~,index]=min(diff_ran);
    % the color of arrow
    colorarrow=colorlist(index,:);
    q =quiver(c(p_i,1),c(p_i,2),v(p_i,1),v(p_i,2),'MaxHeadSize',100);
    q.LineWidth=1;
    q.Color=colorarrow;
    hold on
    s=scatter(c(p_i,1),c(p_i,2),50,'filled');
    s.MarkerFaceColor=colorarrow;
    s.MarkerEdgeColor=colorarrow;
    hold on
end
cl=colorbar;
cl.TickLabels=linspace(0,2*pi,11);
cl.Title.String='\alpha';
cl.Title.FontSize=13;

