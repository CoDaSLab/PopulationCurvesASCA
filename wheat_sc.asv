%
% coded by: Jose Camacho
% last modification: 17/Nov/22

% Make sure MEDA-Toolbox v1.4 is on the path

close all
clear all
load wheat

ut = unique(treatment);
for i=1:length(treatment)
    for j=1:length(ut)
        if strcmp(treatment{i},ut(j))
            ytre(i)=j;
        end
    end
end

ut = unique(time);
for i=1:length(time)
    for j=1:length(ut)
        if strcmp(time{i},ut(j))
            ytim(i)=j;
        end
    end
end

ut = unique(trait);
for i=1:length(trait)
    for j=1:length(ut)
        if strcmp(trait{i},ut(j))
            ytra(i)=j;
        end
    end
end

ut = unique(replicate);
for i=1:length(replicate)
    for j=1:length(ut)
        if strcmp(replicate{i},ut(j))
            yrep(i)=j;
        end
    end
end

F=[ytim' ytra' ytre'];

%% ASCA model: Full experimental matrix
[T, parglmo] = parglm(X, F, 'interaction', 1);

T.Source{2} = 'Time-1';
T.Source{3} = 'Trait-2';
T.Source{4} = 'Treatment-3';

disp(T)

M = size(X,2);

Xs.N = size(X,1);
Xs.M = size(X,2);
Xs.k = [std(parglmo.factors{1}.matrix,0,'all')/sqrt(M),...
    std(parglmo.factors{2}.matrix,0,'all')/sqrt(M),...
    std(parglmo.factors{3}.matrix,0,'all')/sqrt(M),...
    std(parglmo.interactions{1}.matrix,0,'all')/sqrt(M),...
    std(parglmo.interactions{2}.matrix,0,'all')/sqrt(M),...
    std(parglmo.interactions{3}.matrix,0,'all')/sqrt(M)];


[PCmean,PCrep,struct] = powercurve(Xs,F,{[1 2],[1,3],[2,3]},2,1000,[],@()1,[],.05,1,200,1,[],[],[]);  
legend('Time-1','Trait-2','Treat-3','1-2','1-3','2-3')
title('Absolute Sample Curves, 5-Replicate')

%% ASCA model: Single replicate experimental matrix
[B, ia, ic] = unique(F,'rows','stable');

Xc = X(ia,:);
Fc = F(ia,:);
M = size(Xc,2);

[T2, parglmo] = parglm(Xc,Fc,'interaction',1);

T2.Source{2} = 'Time-1';
T2.Source{3} = 'Trait-2';
T2.Source{4} = 'Treatment-3';

disp(T2)




Xs.N = size(X,1);
Xs.M = M;
Xs.k = [std(parglmo.factors{1}.matrix,0,'all')/sqrt(M),...
    std(parglmo.factors{2}.matrix,0,'all')/sqrt(M),...
    std(parglmo.factors{3}.matrix,0,'all')/sqrt(M),...
    std(parglmo.interactions{1}.matrix,0,'all')/sqrt(M),...
    std(parglmo.interactions{2}.matrix,0,'all')/sqrt(M),...
    std(parglmo.interactions{3}.matrix,0,'all')/sqrt(M)];


[PCmean,PCrep,struct] = powercurve(Xs,Fc,{[1 2],[1,3],[2,3]},2,1000,[],@()1,[],.05,1,200,1,[],[],[]);  
legend('Time-1','Trait-2','Treat-3','1-2','1-3','2-3')
title('Absolute Sample Curves, 1-Replicate')



