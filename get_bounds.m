function [lb,ub]=get_bounds()
load('iJO1366.mat')
[~,size0]=size(iJO1366.S);
model = convertToIrreversible(iJO1366);
n_r=size(model.S,2);


%Reversibility

rev_match=zeros(size0,2);
for i=1:size0
    if model.match(i)==0
        rev_match(i,:)=0;
    else
        rev_match(i,1)=i;
        rev_match(i,2)=model.match(i);
    end
end   
    
%Upper and lower bound from media and growth conditions
media=xlsread("Copy of Davidi_media.xlsx");
lb=zeros(n_r,31);
ub=zeros(n_r,31);

for i=1:31
    lb(:,i)=model.lb;
    ub(:,i)=model.ub;
    if rev_match(media(i,1),1)~=0
        ub(rev_match(media(i,1),2),i)=10;
        ub(rev_match(media(i,1),1),i)=1000;
    end
    limit=setdiff(media(:,1),media(i,1)); % media(i,1)=uptake
    lb(limit,i)=0;
    ub(limit,i)=0;

    lb(model.c==1,i)=media(i,2)-0.01; %model.c==1 => biomass index
    ub(model.c==1,i)=media(i,2)+0.01; %media(:,2)= growth rate
end

end