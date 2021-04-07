function [milp]=getkmax(Kapp_matrix,V_matrix,homo)
media=xlsread("Copy of Davidi_media.xlsx");
[k_max,ind_max]=max(Kapp_matrix,[],2);

ind=find(k_max~=0);
milp.kmax=k_max(ind);
milp.reac=string(homo.reac(ind));
milp.reacind=homo.reacind(ind);
milp.genes=string(homo.genes(ind));
milp.v=V_matrix(ind,:);
milp.kapp=Kapp_matrix(ind,:);
milp.conditions=media(ind_max(ind));
milp.reac_rev=homo.reac_rev(ind);
milp.kapp=Kapp_matrix(ind,:);
milp.abun=homo.abun(ind,:);
end