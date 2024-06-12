%
% coded by: Jose Camacho
% last modification: 17/Nov/22

% Make sure MEDA-Toolbox v1.4 is on the path

close all
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

%% ASCA model: Single replicate experimental matrix
[B, ia, ic] = unique(F,'rows','stable');

[T2, parglmo] = parglm(X(ia,:),F(ia,:), 'interaction',1);

T2.Source{2} = 'Time-1';
T2.Source{3} = 'Trait-2';
T2.Source{4} = 'Treatment-3';

disp(T2)
