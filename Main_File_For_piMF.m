%%% Main File for pre-processing and calling the piMF function to recover
%%% reaction networks and the extents of reaction.
%%% An example reaction of transesterification of biodiesel reaction (S.C. Burnham et al. / Chemical Engineering Science 63 (2008) 862–873) 
%%% is used to demonstrate the algorithm.  
%%% data used for illustration: biodiesel_data.mat
%%% contains time series data and atomic matrix for the different chemicals
%%% involved in the system as follows:
%%% time: 31 x 1 dimensional matrix of time stamps
%%% C   : 31 x 6 dimensional matrix of concentration for seven species
%%% A   : 4  x 6 dimensional matrix of atomic matrix

clear all
clc
close all

%% Load data
load biodiesel_data.mat;

%% Initial data analysis

%%% size of data (observations (m_C) x number of species (n_C))
[m_C, n_C] = size(C);

%%% Computing Reaction variant form of C
C_RV = (C-ones(m_C,1)*C(1,:));

%%% Maximum number of reactions
temp = rank(A);
R_max = n_C-rank(A);

%%% Observed number of reactions through singular value decomposition
s_obs = svd(C_RV);
[m_s n_s] = size(s_obs);

%%% Number of significant eigenvalues is computed using 95% thresholding
R = 0;
for i = 1:max([R_max m_s])
    s_thres = (sum(s_obs(1:i))/sum(s_obs))*100;
    if s_thres > 95
        R = i;
        break;
    end 
end

%% Set simulation conditions for piMF and run simulations

%%% Value that stoichiometric coefficient can take
S_set = [-1, 0, 1];

% Sparsity: Maximum number of species that can participate
% set K value is used for this dataset
K = 6;

%%% Search diameter for sphere decoding algorithm
dia_sphere_decoding = 0.5;

%%%% Initialization of algorithm
% for user-defined instants of random initialization recovered extents of reaction
% and stoichiometry 

no_rand_int = 100; % # of random initialization

if ~isempty(gcp('nocreate'))
    delete(gcp('nocreate'));
end

numWorkers = 10; % Specify the number of workers you want
parpool('local', numWorkers);

parfor ig = 1:no_rand_int
W_initial = sort(abs(randn(m_C,R)), 'ascend');

%%%% Calling piMF 
[X_out{ig},N_out{ig},fvalopt(ig)] = piMF(C_RV,S_set,R,A,K,dia_sphere_decoding,W_initial);

end

delete(gcp('nocreate'));

%% Post analysis of piMF results

%%%% finding unique sets of X_out and N_out (removing repetitive results) 

normalizedMatrices = cellfun(@(M) sortrows(M), N_out, 'UniformOutput', false); % Normalize rows of each matrix

serializedMatrices = cellfun(@(M) mat2str(M), normalizedMatrices, 'UniformOutput', false); % Serialize each matrix

[~, uniqueIdx] = unique(serializedMatrices); % Find unique serialized matrices

N_uni_sets = N_out(uniqueIdx); % Retrieve unique sets of stoichiometric matrix

%%%% Computing the errors from piMF 
for i = 1:length(N_uni_sets)
Obj = C_RV - X_out{uniqueIdx(i)}*N_out{uniqueIdx(i)};
fval_pimf(i) = (norm(Obj,'fro'))^2;
end

%%%% find physically feasible stoichiometric matrix and its extents

% constraints used here is species A should be reactant or nonparticpitant
% in all reaction

idx_infeas_N = find(cellfun(@(M) any(M(:,1)>0), N_uni_sets));
idx_feas_N = setdiff(1:length(N_uni_sets),idx_infeas_N); % physically feasible N

N_best = N_uni_sets{idx_feas_N}; 
X_best = X_out{uniqueIdx(idx_feas_N)};

%%%% Arrange the rows of feasible stoichiometric matrix for better convenience
% R1 - A and B reactants, R2 - B and C reactants, R3 - B and E reactants

idx_R1 = find(N_best(:,1)<0);
idx_R2 = find((N_best(:,2)<0)&(N_best(:,3)<0));
idx_R3 = find((N_best(:,2)<0)&(N_best(:,5)<0));

temp_N = [N_best(idx_R1,:);N_best(idx_R2,:);N_best(idx_R3,:)];
temp_X = [X_best(:,idx_R1) X_best(:,idx_R2) X_best(:,idx_R3)];

X_best = temp_X;
N_best = temp_N;

%%% plot for time profile of computed extent
plot(time,X_best,'o')
hold on


%%%% save the pimf results as mat file

save('piMF_res.mat','N_best','X_best','C','time','fval_pimf','N_uni_sets','idx_feas_N');