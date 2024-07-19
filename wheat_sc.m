%% EXAMPLE of Absolute Population Curves with the Wheat data set in 
% "Population Power Curves in ASCA with Permutation Testing". Submitted to 
% the Special Issue in honor of Prof. Age Smilde in his retirement.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
%       Michael Sorochan Armstorng (mdarmstr@ugr.es)
% last modification: 19/Jul/2024
%
% Copyright (C) 2024  University of Granada, Granada
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Select data with all cultivars but CM

close all
clear all
load wheat

ut = unique(trait);
for i=1:length(trait)
    for j=1:length(ut)
        if strcmp(trait{i},ut(j))
            ytra(i)=j;
        end
    end
end

ut = unique(treatment);
for i=1:length(treatment)
    for j=1:length(ut)
        if strcmp(treatment{i},ut(j))
            ytre(i)=j;
        end
    end
end
ytre = ytre(find(ytra~=3))';

ut = unique(time);
for i=1:length(time)
    for j=1:length(ut)
        if strcmp(time{i},ut(j))
            ytim(i)=j;
        end
    end
end
ytim = ytim(find(ytra~=3))';

ut = unique(replicate);
for i=1:length(replicate)
    for j=1:length(ut)
        if strcmp(replicate{i},ut(j))
            yrep(i)=j+5*(ytra(i)-1);
        end
    end
end
yrep = yrep(find(ytra~=3))';

X = X(find(ytra~=3),:);
F = [ytre ytim yrep];


%% ASCA model: Single replicate experimental matrix

Xc = X(find(yrep==1 | yrep==6 | yrep==16),:);
Fc = F(find(yrep==1 | yrep==6 | yrep==16),:);

[Ts, parglmo] = parglm(Xc, Fc, 'interaction', 1, [], 2, [], [], [], [1 3]);
 
Ts.Source{2} = 'Treatment-A';
Ts.Source{3} = 'Time-B';
Ts.Source{4} = 'Individual-C(A)';
Ts.Source{5} = 'Int-AB';

disp(Ts)

MSA = Ts{2,5};
MSB = Ts{3,5};
MSCA = Ts{4,5};
MSAB = Ts{5,5};
MSE = Ts{6,5};


%% APC from MSE estimates

Xs.N = size(Xc,1);
Xs.M = size(Xc,2);
Xs.k = [sqrt(max(MSA+MSE-MSAB-MSCA,0)/(5*5)),...
        sqrt(max(MSB-MSAB,0)/(2*5)),...
        sqrt(max(MSCA-MSE,0)/5),...
        sqrt(max(MSAB-MSE,0)/5),...
        sqrt(MSE)];
    

[PCmean,PCrep,struct] = powercurve(Xs,Fc, 'interaction', 2, [], [], @()1, 1:15, [], 1, [], 2, [], [], [], [1,3], 3);  

legend('Treatment-A','Time-B','Individual-C(A)','Int-AB') 
saveas(gcf,'./Figures/APCW'); saveas(gcf,'./Figures/APCW.eps','epsc'); 


%% ASCA model: Full experimental matrix

[T, parglmo] = parglm(X, F, 'interaction', 1, [], 2, [], [], [], [1 3]);
 
T.Source{2} = 'Treatment-A';
T.Source{3} = 'Time-B';
T.Source{4} = 'Individual-C(A)';
T.Source{5} = 'Int-AB';

disp(T)