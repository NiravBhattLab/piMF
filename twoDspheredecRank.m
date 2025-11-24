  function [Sol_Matrix,nvisit,All_solutions] = twoDspheredecRank(A,X_r,D,nz,Dia,S)
  %%% Function provides solution to constrained integer least-squares
  %%% problem. 
  %%% Inputs %%%
  
[r,row]=size(X_r);
[r,col]=size(A);
col = size(D,2);
%Solution to Ax=0, A*N'=0
[ax_0,nvisit] = get_vec_hnf(A,S,nz);
%ax_0(1)
%G= containers.Map('KeyType','int32','ValueType','any');
G= containers.Map(int32(1),[0 0]);
remove(G,1);


%As many rows in N, that many times tree has to be duplicated
for i=1:row
    G(i)=[ax_0;containers.Map()];
end
AX_0=ax_0;

%start from first column and solve using sphere decoder
allkeys=cell2mat(ax_0.keys());  %Keep a Copy of the solution
sol_mat=cell2mat(ax_0.values()'); %The complete solution matrix of all the rows
min_keys={};                      %Location of solution with min norm    
dia=Dia;
First_col_soln=zeros(row,1);

while length(min_keys)<1 % Run the loop until atleast one solution is found
    %sp_set= containers.Map('KeyType','int32','ValueType','any');
    %solu_set= containers.Map('KeyType','int32','ValueType','any');
    sp_set= containers.Map(int32(1),[0 0]);
    solu_set= containers.Map(int32(1),[0 0]);
    remove(sp_set,1);
    remove(solu_set,1);
    %disp(['Sphere Diameter - ',num2str(dia)]);
    for j=1:col
        %col=2;
%         disp(['Processing Column number - ',num2str(j)]);
        if isempty(solu_set)                                %if solution set is empty
            for i=1:row                                     
                sol_mat_g=cell2mat(G(i).values()');         
                if ~isempty(sol_mat_g)
                    sp_set(row-i+1)=unique(sol_mat_g(:,j)');%Here the set is stored in reverse order since spheredec starts from bottom
                end
            end
            %dia=0.5;
            
            Sp_Vector  = oneDspheredecRank(sp_set,D(:,j),X_r,dia); %Find the first column of matrix using sphere decoder
            Vector=[];
            for i=1:size(Sp_Vector,2)
                [,indx]=ismember(Sp_Vector(:,i)',First_col_soln','rows'); %If the column is already present dont consider it
                if ~indx
                    First_col_soln=[First_col_soln Sp_Vector(:,i)];     %else add it to solution and the list
                    Vector=[Vector Sp_Vector(:,i)];                     %Vector contains first column on which other columns will  be built
                end
            end
            if isempty(Vector)
%                 disp(['No new solutions. Number of solns in prev round - ',num2str(size(First_col_soln,2))])
                break                                                   %If no new solution increment diameter and try again
            end
            %solu_temp=containers.Map('KeyType','int32','ValueType','any');
            solu_temp=containers.Map(int32(1),[1 1]);
            remove(solu_temp,1);
            n=size(Vector,2);
            for t_j=1:n
                solu_temp(t_j)=Vector(:,t_j);
            end
            %solu_set= containers.Map('KeyType','int32','ValueType','any');
            solu_set=containers.Map(int32(1),[1 1]);
            remove(solu_set,1);
            solu_set=[solu_temp;containers.Map()];
        else
            n_sol=1;
            %solu_temp=containers.Map('KeyType','int32','ValueType','any');
            solu_temp=containers.Map(int32(1),[1 1]);
            remove(solu_temp,1);
            for i=1:length(solu_set.keys())
                Prev_sol = solu_set(i);
                %temp_G= containers.Map('KeyType','int32','ValueType','any');
                temp_G=containers.Map(int32(1),[1 1]);
                remove(temp_G,1);
                for t_i=1:row
                    temp_G(t_i)=[G(t_i);containers.Map()];
                end
                % For all prev solutions delete the rows which dont contain
                % required bit at prev columns
                for t_col=1:size(Prev_sol,2)
                    for v_r=1:row
                        bit=Prev_sol(v_r,t_col);
                        loc=t_col;
                        del_keys={};
                        t_n=1;
                        %temp= containers.Map('KeyType','int32','ValueType','any');
                        temp=containers.Map(int32(1),[1 1]);
                        remove(temp,1);
                        temp=[temp_G(v_r);containers.Map()];
                        temp_keys=cell2mat(temp.keys());
                        for t_i=1:length(temp_keys)
                            z=temp(temp_keys(t_i));
                            loc;
                            if z(loc)~=bit
                                del_keys{t_n}=temp_keys(t_i);
                                t_n=t_n+1;
                            end
                            
                        end
                        del_keys;
                        remove(temp_G(v_r),del_keys);
                        %del_keys={};
                    end
                end
                
                %for i=1:G.length()
                %    cell2mat(G(i).values()')
                %end
                %%% After deleting for next column find the set
                length_set=1;
                for t_i=1:row
                    sol_mat_g=cell2mat(temp_G(t_i).values()');
                    if ~isempty(sol_mat_g)
                        sp_set(row-t_i+1)=unique(sol_mat_g(:,j)');%Here the set is stored in reverse order since spheredec starts from bottom
                        length_set=length_set*length(sp_set(row-t_i+1));
                    end
                end
                %%% If set contains more than one elemnt run sphere dec
                if length_set>1
                    %disp(['Processing Column number - ',num2str(j)])
                    Vector  = oneDspheredecRank( sp_set,D(:,j),X_r,dia);
                else
                %%%% Else output set in reverse dir
                    Vector=flipud(cell2mat(sp_set.values()'));
                end
                
                n=size(Vector,2);
                for t_j=1:n
                    solu_temp(n_sol)=[solu_set(i) Vector(:,t_j)];
                    solu_temp.keys();
                    n_sol=n_sol+1;
                end
            end
            %solu_set= containers.Map('KeyType','int32','ValueType','any');
            solu_set=containers.Map(int32(1),[1 1]);
            remove(solu_set,1);
            solu_set=[solu_temp;containers.Map()];
            
            
        end
    
    end
   
  % All_solutions=temp_G;
    

    % Among solutions find the one with lowest norm and full rank
    %solu_set(1)
    min_val = 1e6;
    minind=1;
    min_keys={};
    %solution=containers.Map('KeyType','int32','ValueType','any');
    solution=containers.Map(int32(1),[1 1]);
    remove(solution,1);
    solution=[solu_set;containers.Map()];
    for i=1:length(solu_set.keys())
        if rank(solu_set(i))~=row
            %i,solu_set(i)
            remove(solu_set,i);
            continue
        end
        %i,solu_set(i),
        f_norm=norm(D-X_r*solu_set(i));
        if f_norm<min_val
            %minind=i;
            min_val=f_norm;
            min_keys{minind}=i;
        else
            if abs(f_norm-min_val)<1e-6
                minind=minind+1;
                min_keys{minind}=i;
                
            end
        end
    end
    dia=dia+1;
end
Sol_Matrix = solu_set(min_keys{1});
%fprintf('Total Nodes visited %d\n',nvisit);
%min_keys
%solu_set(minind)