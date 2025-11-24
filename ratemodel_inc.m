% Model candidates for incremental identification

function rate = ratemodel_inc(C,i,j,par)

% Possible model candidates

if(i==1) && (j==1)
    rate = par.*ones(length(C),1); %r1,1
elseif(i==1) && (j==2)
    rate = par.*C(:,1); %r1,2
elseif(i==1) && (j==3)
    rate = par.*C(:,2); %r1,3
elseif(i==1) && (j==4)
    rate = par.*(C(:,1).^2); %r1,4
elseif(i==1) && (j==5)
    rate = par.*(C(:,2).^2); %r1,5
elseif(i==1) && (j==6)
    rate = par.*C(:,1).*C(:,2); %r1,6
elseif(i==1) && (j==7)
    rate = par(1).*C(:,1).*C(:,2)-par(2).*C(:,3).*C(:,4); %r1,7    
    
elseif(i==2) && (j==1)
    rate = par.*ones(length(C),1); %r2,1
elseif(i==2) && (j==2)
    rate = par.*C(:,2); %r2,2
elseif(i==2) && (j==3)
    rate = par.*C(:,3); %r2,3
elseif(i==2) && (j==4)
    rate = par.*(C(:,2).^2); %r2,4
elseif(i==2) && (j==5)
    rate = par.*(C(:,3).^2); %r2,5
elseif(i==2) && (j==6)
    rate = par.*C(:,3).*C(:,2); %r2,6
elseif(i==2) && (j==7)
    rate = par(1).*C(:,3).*C(:,2)-par(2).*C(:,4).*C(:,5); %r2,7     
    
elseif(i==3) && (j==1)
    rate = par.*ones(length(C),1); %r3,1
elseif(i==3) && (j==2)
    rate = par.*C(:,2); %r3,2
elseif(i==3) && (j==3)
    rate = par.*C(:,5); %r3,3
elseif(i==3) && (j==4)
    rate = par.*(C(:,2).^2); %r3,4
elseif(i==3) && (j==5)
    rate = par.*(C(:,5).^2); %r3,5
elseif(i==3) && (j==6)
    rate = par.*C(:,2).*C(:,5); %r3,6
elseif(i==3) && (j==7)
    rate = par(1).*C(:,2).*C(:,5)-par(2).*C(:,4).*C(:,6); %r3,7     
    
end

