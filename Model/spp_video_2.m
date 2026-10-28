clear variables
clc

% the number of individuals
N=50;
% the weighting term for individuals
omega=0.4;
% the preferred direction of informed individuals
g=[2 3];
g=g./sqrt(sum(g.^2));
% the length of time step
delta_t=0.1;
% the time of evolution
T=1000;
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

rng(5,'twister')
% the initial condition: random positions and orientations
s_lat=3;
c_x=s_lat/2*(2*rand(N,1)-1);c_y=s_lat/2*(2*rand(N,1)-1);
% the unit vector
ran_angle=2*pi*rand(N,1);
v_x=cos(ran_angle);v_y=sin(ran_angle);
c=[c_x c_y]; v=[v_x v_y];

% the size of the uncertainty range within the two-component group, the
% control parameter is denoted by p, the proportion of informed individuals
p=0.2;
% the uncertainty range
ran=2*pi*ones(N,1);
ran(1:round(N*p))=0;

% the center of mass
cent=zeros(Nt,2);
% the deviation of an individual's preferred direction from the
% accurate direction
dev=zeros(Nt*N,1);

% the video file
spp_v='spp_video_2.avi';
Vid_Wri_2=VideoWriter(spp_v);
Vid_Wri_2.FrameRate=10;
open(Vid_Wri_2);

f=figure('color','white');


for t=1:Nt
    
    new_v=zeros(N,2); % all individuals require update simultaneously!
    if t>550 && t<650
        q =quiver(c(1:round(N*p),1),c(1:round(N*p),2),v(1:round(N*p),1),v(1:round(N*p),2),0.5);
        q.LineWidth=1;
        q.Color=[1 0 0];
        hold on
        ss=scatter(c(1:round(N*p),1),c(1:round(N*p),2),30,'filled');
        ss.MarkerFaceColor=[1 0 0];
        ss.MarkerEdgeColor=[1 0 0];
        hold on
        q =quiver(c(round(N*p)+1:end,1),c(round(N*p)+1:end,2),v(round(N*p)+1:end,1),v(round(N*p)+1:end,2),0.5);
        q.LineWidth=1;
        q.Color=[0 0 1];
        hold on
        ss=scatter(c(round(N*p)+1:end,1),c(round(N*p)+1:end,2),30,'filled');
        ss.MarkerFaceColor=[0 0 1];
        ss.MarkerEdgeColor=[0 0 1];
        xlabel('{\it x}','fontsize',15,'FontName','TimesNewRoman')
        ylabel('{\it y}','fontsize',15,'FontName','TimesNewRoman')
        hold off
    
        frame=getframe(gcf);
        writeVideo(Vid_Wri_2,frame);
        pause(0.05)
    end
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

close(Vid_Wri_2)