function v_hat=topological_int(v,c,N,k,i,alpha)

        % this function calculates the social interactions acting on an individual
        % by its near neighbours. The key point is that interacting neighbours is
        % determined by topological rules, namely, k-nearest neighbours and evenly
        % distributed in the two sides with respect to the direction of motion of
        % the focal individual
        
        % the input:
        % v and c, positions and velocities of all individuals
        % N, the group size
        % k, the k-nearest neighbours
        % i, the focal individual
        % alpha, the repulsion range
        
        % the output:
        % v_hat, the social interaction acting on the focal individual
        
        % define two summation vectors to calculate the velocity according to the
        % repulsion and attraction rules, respectively
        v_r=[0 0];
        v_a=v(i,:);
        
        % calculating the repulsion interaction
        for j=1:N
            dis=sqrt(sum((c(i,:)-c(j,:)).^2));
            if dis>0 && dis<=alpha
               v_r=v_r-(c(j,:)-c(i,:))./sqrt(sum((c(j,:)-c(i,:)).^2));
            end
        end
        
        % determine the social interaction acting on the focal individual 
        if sum(v_r.^2)>0
            v_hat=v_r./sqrt(sum(v_r.^2));
        else
            % an angle threshold that limits the spatial positions of candidate 
            % interacting neighbours
            ang_thre=pi;

            % the angle between the direction of motion of the focal individual
            % and the relative position of all other individuals with respect
            % to the focal individual
            ang=acos(sum(v(i,:).*(c-c(i,:)),2)./(sqrt(sum(v(i,:).^2))....
                *sqrt(sum((c-c(i,:)).^2,2))));
            % the vector product between the above two vectors
            vec_pro=cross([ones(N,2).*v(i,:),zeros(N,1)],[c-c(i,:),zeros(N,1)]);
            % define two vector to record the neighbours of the focal
            % individual in the right and left sides, respectively
            right_nei=zeros(N,1); right_s=0;
            left_nei=zeros(N,1); left_s=0;
            for j=1:N
                if ang(j)<ang_thre
                    if vec_pro(j,3)<0
                        right_s=right_s+1;
                        right_nei(right_s)=j;
                    elseif vec_pro(j,3)>0
                        left_s=left_s+1;
                        left_nei(left_s)=j;
                    end
                end
            end
            right_nei=right_nei(1:right_s);
            left_nei=left_nei(1:left_s);
            % calculating the distance between the focal individuals and the
            % candidate neighbours
            right_nei_dis=sqrt(sum((c(i,:)-c(right_nei,:)).^2,2));
            left_nei_dis=sqrt(sum((c(i,:)-c(left_nei,:)).^2,2));
            % there may be the situation that the candidate neighbours is not enough to
            % reach k/2, we therefore need to determine the number of interacting 
            % neighbours in left and right sides, respectively
            if length(left_nei)<k/2
                N_left=length(left_nei);
            else
                N_left=k/2;
            end
            if length(right_nei)<k/2
                N_right=length(right_nei);
            else
                N_right=k/2;
            end
            % ranking the distance in order from small to large
            [~,Mr]=sort(right_nei_dis);
            [~,Ml]=sort(left_nei_dis);
            for kr=1:N_right
                r_index=right_nei(Mr(kr));
                v_a=v_a+(c(r_index,:)-c(i,:))./sqrt(sum((c(r_index,:)-c(i,:)).^2))+...
                     v(r_index,:)./sqrt(sum(v(r_index,:).^2));
            end
            for kl=1:N_left
                l_index=left_nei(Ml(kl));
                v_a=v_a+(c(l_index,:)-c(i,:))./sqrt(sum((c(l_index,:)-c(i,:)).^2))+...
                     v(l_index,:)./sqrt(sum(v(l_index,:).^2));
            end

            v_hat=v_a./sqrt(sum(v_a.^2));

        end

end
