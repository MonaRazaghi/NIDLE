function [Kapp_matrix,V_matrix]=getkapp_MILP(R_index,R_ab_matrix,sol_flux,sol_y)

Kapp_matrix=zeros(length(R_index),31);
V_matrix=zeros(length(R_index),31);


for cond=1:31
    y=sol_y(:,cond);
    V=sol_flux(:,cond);  
    active=R_index(y>0.999 & y<1.001);

    for i=1:length(R_index)
        if find(R_index(i)==active)
            Kapp_matrix(i,cond)=V(R_index(i))/3600/R_ab_matrix(i,cond);
            V_matrix(i,cond)=V(R_index(i));
            % if the abundance is not NAN, and flux is not zero, it is a
            % used enzyme
            
        end
    end
end
end