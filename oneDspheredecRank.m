function [ Vector ] = oneDspheredecRank( set,y,H,D )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[N K]=size(H);
[Q R]=qr(H); %QR
yhat = Q'*y; %multiply by Q'
%i=N
Vec_prev =[];
nz_col = find(R(end,:)~=0);
if ~isempty(nz_col)
    M=length(nz_col);
    Z=genAllPossibleVectors_multiD(set,M,1);
    disp(['Number of candidate vectors of dimension - ',num2str(M),' are ',num2str(size(Z,2))]);
end
while isempty(Vec_prev) %increase diameter if no candidate found
    prev_row_len=0;
    Vec_prev =[];
    row_idx=max(find(sum(abs(R'))));%Find index of first non-zero row
    for i = row_idx:-1:1 %start from last row
        %disp(['Processing row number - ',num2str(i)])
        nz_col = min(find(R(i,:))); %starting index of the first non-zero element in ith row of R
        M=length(R(i,nz_col:end));  %number of non-zero elements in ith row of R
        coeff= R(i,nz_col:end); %Get the coefficients
        
        Z=genAllPossibleVectors_multiD(set,M-prev_row_len,row_idx-i+1); % Gen all possible vectors starting from dimension = row_idx-i+1 upto M-prev_row_len dimensions
        % Since the row processing is from last and the set passed has last
        % row in first hence the reversal
        Vec_row=[];
        
        if ~isempty(Vec_prev)
            
            for j=1:size(Vec_prev,2)
                for k=1:size(Z,2)
                    one_vec=[Z(:,k);Vec_prev(:,j)];
                    radius=yhat(i)-coeff*one_vec;
                    Dia = radius.^2;
                    if Dia <= D
                        Vec_row=[Vec_row one_vec];
                    end
                end
            end
            
        else
            oney = ones(1,size(Z,2));
            y2 = yhat(i)*oney;
            %size(y2)
            %size(coeff*Z)
            Dia = diag((y2 - coeff*Z)'*(y2-coeff*Z));
            diaind = Dia<=D;
            Vec_row = Z(:,diaind);
        end
        
        %Vec_row
        %radius=yhat(i)-coeff*Vec_row %(yhat - h_ix_i)^2
        %Dia = radius.^2
        %index = find(Dia<=D) %choose vectors with dia less than D
        if ~isempty(Vec_row)>0
            Vec_prev=Vec_row;%(:,index)
            prev_row_len=M;
        else
            Vec_prev=[];
            break
        end
        %disp(['Number of vectors retained is ',num2str(length(Vec_prev))]);
    end
    D=D+1;
    %Vec_prev=[]
end
%Vec_row
radius=yhat(i) - coeff*Vec_row;
Dia=radius.^2;
%[minv mind]=min(Dia);
minv=min(Dia);
%minds=find(Dia<=minv+0.01);
%Vec_prev
%M_dist_vec= Vec_prev(:,minds);
M_dist_vec= Vec_prev; %Return all vectors with diameter
Vector = M_dist_vec;
end


function d1 = genAllPossibleVectors_multiD(multiDdict,K,idx)
%usage : d1= genAllPossibleVectors(x1,K)
% Generates all possible vectors of length K containing elements of the set 
% multiDdict  is the dictionary containing set of elements for each
% dimension
% K is the number of dimension
%idx is the starting dimension
% this function returns x1^K
%last=length(multiDdict);
x1=multiDdict(idx);
d1 = x1;
for k=idx+1:idx+K-1
    temp = d1;
    clear d1;
    %d1=[]
    x1=multiDdict(k);
    for k1=1:length(temp)
       for k2=1:length(x1)
           index=(k1-1)*length(x1) + k2;
           %index,k,k1,k2;
           d1(:,(k1-1)*length(x1) + k2) = [x1(k2);temp(:,k1)];
       end;
    end;
end;
end