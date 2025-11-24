function [X_out,N_out,fvalopt]=piMF(C_RV,S_set,R,A,K,Dia,W_in)
%%%% piMF provides the extents of reactions (X) and reaction network (N)
%%%%% Inputs %%%%%
% C_RV: n X m dimensional reaction variants concentration matrix
% S_set: Set of integer the stoichiometric matrix can take as a vector from
% the lowest to the largest value
% R : Number of reactions in the system or rank of reaction network (N)
% K:  Sparsity: Maximum number of species that can participate
% Dia: Search diameter for sphere decoding algorithm
% W_in: Initialization of the extents of reaction 
%%%%% Outputs %%%%
% X_out: n x R dimensional matrix of the extents of reaction
% N_out : R x m dimensional matrix of the reaction networks
%%%%%%%%%%% 
tol = 10^-6;
MaxIter = 10^3; %max iter of While loop

%%% fmincon options
MaxFunEvals = 100;
MaxIterFminCon = 100;
TolFun = 10^-6;
%%% Initialization and Computation of necessary matrices
    [n m]=size(C_RV);
    r1=[1 -1 zeros(1,n-2)];
    c1=[1;zeros(n-2,1)];
    T=toeplitz(c1,r1);
    C2=T(1:n-1,:);
    I=eye(R);
    Aw=kron(I,C2);
    bw=tol*ones(R*(n-1),1);
    Xold =ones(n,R);
    Nold = ones(R,m);
    X = W_in;
    N = randi([S_set(1) S_set(end)],R,m);
%%%

options = optimset('display','iter','MaxFunEvals',MaxFunEvals,'MaxIter',MaxIterFminCon,'TolFun',TolFun);
k = 0;

while  k < MaxIter  && (sum(sum((Nold - N).^2)) > tol  || sum(sum((Xold - X).^2)) > tol)
%while  (k < MaxIter &&  all(N(:,1) <= 0) ||  all(N(:,1) ~= 0) &&  all(N(:,2) <= 0) || all(N(:,2) ~= 0)) && (sum(sum((Nold - N).^2)) > tol  || sum(sum((Xold - X).^2)) > tol)
    ht=sum(sum((Nold - N).^2));
    wt=sum(sum((Xold - X).^2));
    Xold = X;
    Nold = N;
%%% Estimation of N given X  using the proposed algorithm in Ayyar and Bhatt, 2019  
    N=twoDspheredecRank(A,X,C_RV,K,Dia,S_set);
%%% Estimation of X given N
    P=kron(N',eye(n));
    Qw=+2*P'*P;
    fw = - 2*C_RV*N';
    LBw = 1*10^-8*ones(size(X));
    UBw = 500*ones(size(N));
%%% Quadratic programming with constraints
    [w, fval]= quadprog(Qw,fw(:),Aw,bw,[],[],LBw,UBw,X(:),options);
    X = vec2mat(w,size(X,1))';%%
    k = k+1;
end
Nold
N
Xold
X

 ht=sum(sum((Nold - N).^2));
 wt=sum(sum((Xold - X).^2));
X_out = X;
N_out = N;
fvalopt = sum(sum(C_RV- X*N).^2);
end