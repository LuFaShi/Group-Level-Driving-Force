function vel_rot=rotating(vel,vec_plane,theta)

% the input: 
% 1.vel, the velocity vector
% 2.vec_plane, a vector to determine the plane of rotating with the velocity vector itself 
% 3.theta, the angle of rotating

% the output:
% vel_rot, the velocity vector after rotating

% the vector product of the velocity vector and the vec_plane
vec_pro_i=[vel(2)*vec_plane(3)-vel(3)*vec_plane(2) vel(3)*vec_plane(1)-vel(1)*vec_plane(3) ...
    vel(1)*vec_plane(2)-vel(2)*vec_plane(1)];
% normalization
vec_pro_i=vec_pro_i./sqrt(sum(vec_pro_i.^2));
% the vector product of the vec_pro_i and the velocity vector
vec_pro_j=[vec_pro_i(2)*vel(3)-vec_pro_i(3)*vel(2) vec_pro_i(3)*vel(1)-vec_pro_i(1)*vel(3) ...
    vec_pro_i(1)*vel(2)-vec_pro_i(2)*vel(1)];
% the velocity after rotating in a new coordinate system
vel_rot_=[cos(theta) sin(theta) 0];

% transform the vec_rot_ to the original coordinate system
vel_rot=vel_rot_(1)*vel+vel_rot_(2)*vec_pro_j+vel_rot_(3)*vec_pro_i;

end








