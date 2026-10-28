% the control parameter for the four types of groups

% the homogeneous group
b_sets=[0:0.4:5.6,5.9:0.1:6.2];
% the uniform group
r_sets=0:0.4:2*pi;
% the powerlaw group
k_sets=[-0.9:0.1:5.9,6:1:30];
% the two-component group
p_sets=0:0.1:1.0;

% their theoretical mean of the theta-distribution
hom_mean=b_sets./4;
uni_mean=r_sets./8;
pow_mean=(pi/2).*((k_sets+1)./(k_sets+2));
two_com_mean=pi/2.*(1-p_sets);

% their theoretical variance
hom_vari=(b_sets.^2)./48;
uni_vari=7/576.*(r_sets.^2);
pow_vari=(((k_sets+1)./(3*(k_sets+3)))-(((k_sets+1).^2)./(4.*((k_sets+2).^2)))).*pi^2;
two_vari=((-1/4).*(p_sets.^2)+(1/6).*p_sets+1/12).*(pi^2);

% their theoretical two-order approximations of the theta-distribution
hom_two_appro=1-((b_sets.^2)./24);
uni_two_appro=1-((r_sets.^2)./72);
pow_two_appro=1-(pi^2/6).*((k_sets+1)./(k_sets+3));
two_two_appro=1+(pi^2/6).*(p_sets-1);

% their theoretical three-order approximations
hom_three_appro=1-((b_sets.^2)./24)+((b_sets.^4)./1920);
uni_three_appro=1-((r_sets.^2)./72)+((r_sets.^4)./9600);
pow_three_appro=1-(pi^2/6).*((k_sets+1)./(k_sets+3))+(pi^4/120).*((k_sets+1)./(k_sets+5));
two_three_appro=1+(pi^2/6).*(p_sets-1)+(pi^4/120).*(1-p_sets);

% their theoretical four-order approximations
hom_four_appro=1-((b_sets.^2)./24)+((b_sets.^4)./1920)-((b_sets.^6)./322560);
uni_four_appro=1-((r_sets.^2)./72)+((r_sets.^4)./9600)-((r_sets.^6)./2257920);
pow_four_appro=1-(pi^2/6).*((k_sets+1)./(k_sets+3))+...
    (pi^4/120).*((k_sets+1)./(k_sets+5))-(pi^6/5040).*((k_sets+1)./(k_sets+7));
two_four_appro=1+(pi^2/6).*(p_sets-1)+(pi^4/120).*(1-p_sets)+(pi^6/5040).*(p_sets-1);

figure(1)
plot(hom_mean,hom_two_appro);
hold on
plot(uni_mean,uni_two_appro);
plot(pow_mean,pow_two_appro);
plot(two_com_mean,two_two_appro);

figure(2)
plot(hom_mean,hom_three_appro);
hold on
plot(uni_mean,uni_three_appro);
plot(pow_mean,pow_three_appro);
plot(two_com_mean,two_three_appro);

figure(3)
plot(hom_mean,hom_four_appro);
hold on
plot(uni_mean,uni_four_appro);
plot(pow_mean,pow_four_appro);
plot(two_com_mean,two_four_appro);