function [new]=MILP(model_orig,homo,lb,ub)


%split the reversible reactions into 2 irreversible reactions
model = convertToIrreversible(model_orig); 
[n_m,n_r]=size(model.S);

new.y=zeros(length(homo.reacind),31);
new.flux=zeros(n_r,31);
new.active=zeros(length(homo.reacind),31);


for cond=1:size(lb,2)
    
    %Assign y variables to reactions with not nan abundance
    ind_notnan_inv=homo.reacind(~isnan(homo.abun(:,cond)));
    ind_notnan_inh=find(~isnan(homo.abun(:,cond)));
    ind_f=find(contains(homo.reac(ind_notnan_inh),'_f'));
    ind_b=find(contains(homo.reac(ind_notnan_inh),'_b'));

    
    n_y=length(ind_notnan_inh); 
    y0=n_r+1; %index where the y variable starts


    %%Maximizing sum of y

    %a) Inequality matrix
    ep=0.0001;
    %Only Y variable constraints only for reactions associated with genes of known abundance
    %Two inequalities of length(abun_genes) as rows, and n_r+n_y as columns
    %One inequality of length of reversible abundance reactions as rows
    I=eye(n_r);
    I_abun=I(ind_notnan_inv,:);
    I_y=eye(n_y);
    

    Vmin=lb(:,cond);
    Vmax=ub(:,cond);

    A1=[-I_abun, -(Vmin(ind_notnan_inv)-ep).*I_y;
        I_abun, (Vmin(ind_notnan_inv)-ep-Vmax(ind_notnan_inv)).*I_y];

    A2=zeros(length(ind_f),n_y);
    for i=1:length(ind_f)
        A2(i,ind_f(i))=1;
        A2(i,ind_b(i))=1;
    end

    A=[A1;[zeros(length(ind_f),n_r),A2]];

    b=[-Vmin(ind_notnan_inv);ep*ones(n_y,1);ones(length(ind_f),1)];

    Aeq=[model.S,zeros(n_m,n_y)];
    beq=zeros(n_m,1);


    f=[zeros(1,n_r),-ones(1,n_y)];
    int=n_r+1:n_r+n_y;

    lower=[lb(:,cond);zeros(n_y,1)];
    upper=[ub(:,cond);ones(n_y,1)];

    sol = intlinprog(f,int,A,b,Aeq,beq,lower,upper);

    y=sol(y0:end);
    active=ind_notnan_inh(y<1.001 & y>0.999);
    
    
    %Minimizing sum of fluxes
    Aeq2=[model.S, zeros(n_m,n_y);
        zeros(1,n_r), ones(1,n_y)];
    beq2=[zeros(n_m,1);sum(y)];
    
    
    f2=zeros(1,n_r+n_y);
    f2(1:n_r)=1; %minimize the total sum of fluxes
    sol_new=intlinprog(f2,int,A,b,Aeq2,beq2,lower,upper);
    %for first condition: sum(sol_new(1:n_r))=1.2432e+04
    y_new=sol_new(n_r+1:end);
    active2=find(y_new<1.001 & y_new>0.999);
    
    new.flux(:,cond)=sol_new(1:n_r);
    new.y(ind_notnan_inh,cond)=y_new;
    new.active(1:length(active2),cond)=ind_notnan_inh(active2);
     
end
end
