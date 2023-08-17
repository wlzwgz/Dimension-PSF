% Clear the workspace
clear all

% Start a timer
tic

% Read the Data
[numdata,text,alldata] = xlsread('Dissimilarities.xlsx');
FungalPat = nan(18,18);
RootPat   = nan(18,18);
Oomycete  = nan(18,18);
for i = 1:length(numdata(:,2))
    FungalPat(numdata(i,2),numdata(i,3)) = numdata(i,8);
    FungalPat(numdata(i,3),numdata(i,2)) = numdata(i,8);
    RootPat(numdata(i,2),numdata(i,3))   = numdata(i,9);
    RootPat(numdata(i,3),numdata(i,2))   = numdata(i,9);
    Oomycete(numdata(i,2),numdata(i,3))  = numdata(i,10);
    Oomycete(numdata(i,3),numdata(i,2))  = numdata(i,10);
end

% Create an interaction matrix using the empirical regressions
A = (1.59-1.52.*FungalPat-2.27.*Oomycete-1.27.*RootPat); % Set off-diagonal entries
A(isnan(A)==1) = -3.4710;                                % Set diagonal entries (dissimilarity = 1)
A = 0.5.*A;                                              % From (symmetrical) feedback to coefficients

% Create a matrix containing all possible species combinations
Combos         = fullfact(2.*ones(1,18))-1;              % Create vectors with all species combinations
[NoComms NoSp] = size(Combos);                           % Save no of communities and max sp. number

% Create an outputmatrix to store the data
OutputMatrix   = nan(NoComms-1,2);                       % Create the output matrix

% Run analysis through all communities
for i = 2:NoComms                           % First entry: only zeros
    % Constraints: community-level feedback, feasibility, local stability
    SpeciesPres    = find(Combos(i,:)==1);
    [~,NoEqSp]     = size(SpeciesPres);
    NoComposites   = sum(SpeciesPres<=6);
    NoGrasses      = sum((SpeciesPres>6).*(SpeciesPres<=12));
    NoLegumes      = sum((SpeciesPres>12));
    NoFamilies     = (NoComposites>0)+(NoGrasses>0)+(NoLegumes>0);
    if NoEqSp > 1
        [Densities, Eigenvalue1, JacobianNumeric, I_s_tot, EigenVPos, Value3] = equilibrium_JJ(A(SpeciesPres,SpeciesPres));
        if I_s_tot < 0 && min(Densities)> 0 && max(Densities)< 1 && real(Eigenvalue1)<0  
            Pairwise = nan(NoEqSp,NoEqSp); 
            PredPSF  = nan(NoEqSp,NoEqSp);
            TempSpeciesMatrixPrevious = A(SpeciesPres,SpeciesPres);
            for q=1:NoEqSp-1
                for qq = q+1:NoEqSp
                    Pairwise(q,qq)=(NoEqSp.*TempSpeciesMatrixPrevious(q,q)+sum(sum(eye(NoEqSp,NoEqSp).*TempSpeciesMatrixPrevious))...
                    -sum(TempSpeciesMatrixPrevious(q,:))-sum(TempSpeciesMatrixPrevious(:,q)))./(NoEqSp-1);
                    PredPSF(q,qq) = Densities(q).*Densities(qq).*Pairwise(q,qq);
                end
            end
        OutputMatrix(i,1)  = NoEqSp;                                          % Species count
        OutputMatrix(i,2)  = NoFamilies;                                      % Number of families
        OutputMatrix(i,3)  = I_s_tot;                                         % Community-level feedback
        OutputMatrix(i,4)  = nanmean(nanmean(Pairwise));                      % Average pair-wise feedback
        OutputMatrix(i,5)  = nansum(nansum(PredPSF));                         % Predicted PSF
        OutputMatrix(i,6)  = -OutputMatrix(i,5).*(1-(1./OutputMatrix(i,1)));  % Pathogen dilution
        OutputMatrix(i,7)  = mean(TempSpeciesMatrixPrevious*Densities');      % W hat
        OutputMatrix(i,8)  = Eigenvalue1;                                     % Dominant Eigenvalue
        OutputMatrix(i,9)  = NoComposites>=1;                                 % Number of composites
        OutputMatrix(i,10) = NoGrasses>=1;                                    % Number of grasses
        OutputMatrix(i,11) = NoLegumes>=1;                                    % Number of legumes
        end
    
    % Constraints: community-level feedback (and, possibly, local stability)
    % Relativize the interaction matrix
    B = A(SpeciesPres,SpeciesPres) - (sum(A(SpeciesPres,SpeciesPres),2)./length(SpeciesPres)) + mean(mean(A(SpeciesPres,SpeciesPres)));
    [Densities, Eigenvalue1, JacobianNumeric, I_s_tot, EigenVPos, Value3] = equilibrium_JJ(B);
    PairwiseB = nan(NoEqSp,NoEqSp); 
    PredPSFB  = nan(NoEqSp,NoEqSp);
    TempSpeciesMatrixPreviousB = B;
    for q=1:NoEqSp-1
        for qq = q+1:NoEqSp
            PairwiseB(q,qq)=(NoEqSp.*TempSpeciesMatrixPreviousB(q,q)+sum(sum(eye(NoEqSp,NoEqSp).*TempSpeciesMatrixPreviousB))...
            -sum(TempSpeciesMatrixPreviousB(q,:))-sum(TempSpeciesMatrixPreviousB(:,q)))./(NoEqSp-1);
            PredPSFB(q,qq) = Densities(q).*Densities(qq).*PairwiseB(q,qq);
        end
    end
    OutputMatrix(i,12) = NoEqSp;                                          % Species count (relativized)
    OutputMatrix(i,13) = NoFamilies;                                      % Number of families (relativized)
    OutputMatrix(i,14) = I_s_tot;                                         % Community-level feedback (relativized) 
    OutputMatrix(i,15) = nanmean(nanmean(PairwiseB));                     % Average pair-wise feedback (relativized)
    OutputMatrix(i,16) = nansum(nansum(PredPSFB));                        % Predicted pair-wise feedback (relativized)
    OutputMatrix(i,17) = -OutputMatrix(i,14).*(1-(1./OutputMatrix(i,1))); % Pathogen dilution (relativized)
    OutputMatrix(i,18) = mean(TempSpeciesMatrixPreviousB*Densities');     % W hat (relativized)
    OutputMatrix(i,19) = Eigenvalue1;                                     % Dominant eigenvalue (relativized)
    OutputMatrix(i,20) = NoComposites>=1;                                 % Number of composites (relativized)
    OutputMatrix(i,21) = NoGrasses>=1;                                    % Number of grasses (relativized)
    OutputMatrix(i,22) = NoLegumes>=1;                                    % Number of legumes (relativized)
    end
end
% Stop the timer
toc

% Optional: save the output as a .mat file
% save OutputData.mat