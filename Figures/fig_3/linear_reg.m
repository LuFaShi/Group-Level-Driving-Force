function [a,b]=linear_reg(x,y)

if size(x)~=size(y)
    x=x';
end

% the mean of x variable
mean_x=mean(x);
% the mean of y variabel
mean_y=mean(y);

% b
b=sum((x-mean_x).*(y-mean_y))/sum((x-mean_x).^2);
% a
a=mean(y)-b*mean(x);

end