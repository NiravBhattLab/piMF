function [G,nnodesvisit] =get_vec_hnf(A,S,nz)
[H U]=Hermite_normal_form_row(A);
%S=[-1,0,1];
%G= containers.Map('KeyType','int32','ValueType','any');
G=containers.Map(int32(1),[1 0]);
remove(G,1);
[rows,cols]=size(A);
prev_row_len=0;
nnodesvisit=0;
for r=rows:-1:1
    row = H(r,:);
    r1=row(min(find(row)):end);
    %r1= row(row~=0);
    if ~isempty(r1)
        if length(r1)==prev_row_len
            temp=containers.Map(int32(1),[1 0]);
            remove(temp,1);
            n=1;
            nodes=G;
            if ~isempty(nodes.keys())
                for c= 1:length(nodes.keys())
                    if nnz(nodes(c))<=nz
                    prev_node = nodes(c);
                    %r1
                    %prev_node
                    val = r1*prev_node';
                       if val==0
                            temp(n)=prev_node;
                            n=n+1;
                       end
                    end
                end
            end
            G=temp;
            
        else    
            for i=1:length(r1)-prev_row_len-1
                
                [G,nnodes]=extend_tree(G,S,nz,r1,'free');
                
                nnodesvisit=nnodesvisit+nnodes;
%                 fprintf('Nodes visited in this stage %d ',nnodes);
%                 fprintf('| Total Nodes visited till now %d\n',nnodesvisit);
%                 fprintf('Nodes Retained for next stage %d\n',G.Count);
                
            end
            
            [G,nnodes]=extend_tree(G,S,nz,r1,'constraint');
            %G(1)
%             fprintf('Nodes visited in this stage %d ',nnodes);
            nnodesvisit=nnodesvisit+nnodes;
%             fprintf('| Total Nodes visited till now %d\n',nnodesvisit);
%             fprintf('Nodes Retained %d\n',G.Count);
            prev_row_len=length(r1);
            %G.Count
        end
    end
    %r,G.Count
end
end

%row=[ 1    -1    -1    -2    -1    -1];

function [temp,nnodes]=extend_tree(G,S,nz,row,type)
nodes = G;
%G.Count
%temp=containers.Map('KeyType','int32','ValueType','any');
temp=containers.Map(int32(1),[1 0]);
remove(temp,1);
n=1;
nnodes=0;
if ~isempty(nodes.keys())
    nnodes = length(nodes.keys())*length(S);
    switch(type)
        case{'free'}
            for c= 1:length(nodes.keys())
                prev_node = nodes(c);
                
                for j=1:length(S)
                    vec =[S(j) prev_node];
                    %vec=[prev_node S(j)];
                    if nnz(vec)<=nz
                        temp(n)=vec;
                        n=n+1;
                        
                    end
                end
            end
        case{'constraint'}
            for c= 1:length(nodes.keys())
                prev_node = nodes(c);
                vec1 =[prev_node];
                val = row(1)*S+row(2:end)*vec1';
                if prod(val)==0
                    sol = S(find(val==0));
                    vec = [sol vec1];
                    if nnz(vec)<=nz && length(sol)>0
                        
                        temp(n)=vec;
                        n=n+1;
                    end
                end
            end
            
            
    end
    
else
    for j=1:length(S)
        temp(n)=S(j);
        n=n+1;
        
    end
end
end


