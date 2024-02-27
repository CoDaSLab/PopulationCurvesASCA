%% EXAMPLES of NULL distributions and ASCA tables in "Power Curves in ASCA with 
% Permutation Testing". Submitted to the Special Issue in honor of Prof. Age
% Smilde.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 26/Feb/2024
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

%% NULL distributions and ASCA table from data simulated with seed = 1. 
% The design matrix F contains a full factorial design with 
% four levels for A, three levels for B and four individuals in each cell 
% of C(A). Other inputs are M = 400, kA = kB = kC(A) = kAB = 0.1*theta 
% and kE = 1-theta, theta = 0, R = 1000, P = 200 and alpha = 0.05.

clear
close all
clc

reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3],1:reps};
theta = 0.5;
k = [theta*[.1,.1,.1,.1], (1-theta)]; 

F = create_design(levels,1);

rng(1);
    
for f=1:2 % main factors A and B
    Xf{f} = randn(length(levels{f}),vars);
end
Xf{3} = randn(length(levels{1})*length(levels{3}),vars);% C(A)
Xi = randn(length(levels{1})*length(levels{2}),vars); % AB
Xe = randn(size(F,1),vars); % E

% normalize the matrices
for f=1:3
    Xf{f} = sqrt(size(Xf{f},1))*Xf{f}/norm(Xf{f},'fro');
end
Xi = sqrt(size(Xi,1))*Xi/norm(Xi,'fro');
Xe = sqrt(size(Xe,1))*Xe/norm(Xe,'fro');

X = [];  
Xf1 = [];
Xf2 = [];
Xf3 = [];
Xi1 = [];
N = size(F,1); 
for i = 1:N     
    Xf1(i,:) = Xf{1}(F(i,1),:);
    Xf2(i,:) = Xf{2}(F(i,2),:);
    Xf3(i,:) = Xf{3}(F(i,1)+(F(i,3)-1)*length(levels{1}),:);
    Xi1(i,:) = Xi(F(i,1)+(F(i,2)-1)*length(levels{1}),:); 
end

X = k(1)*Xf1 + k(2)*Xf2 + k(3)*Xf3 + k(4)*Xi1 + k(5)*Xe; 

T = parglm_dist(X, F, {[1 2]},1,[],2,[],[],[],[1 3]);
legend('A','B','C(A)','AB')
axis([0.8 1.2 0 200])
saveas(gcf,'./Figures/Null3_0'); saveas(gcf,'./Figures/Null3_0.eps','epsc'); 

T.Source(2:5) = {'A','B','C(A)','AB'} 
table2latex(T,'./Figures/Table3_0.tex')

%% NULL distributions and ASCA table from data simulated with seed = 1. 
% The design matrix F contains a full factorial design with 
% four levels for A, three levels for B and four individuals in each cell 
% of C(A). Other inputs are M = 400, kA = kB = kC(A) = kAB = 0.1*theta 
% and kE = 1-theta, theta = 0, R = 1000, P = 200 and alpha = 0.05. The 
% whole experiment is duplicated.

clear
close all
clc

reps = 4;
vars = 400;
levels = {[1,2,3,4],[1,2,3],1:reps};
theta = 0.5;
k = [theta*[.1,.1,.1,.1], (1-theta)]; 

F = create_design(levels,1);
F = [F;F];

rng(1);
    
for f=1:2 % main factors A and B
    Xf{f} = randn(length(levels{f}),vars);
end
Xf{3} = randn(length(levels{1})*length(levels{3}),vars);% C(A)
Xi = randn(length(levels{1})*length(levels{2}),vars); % AB
Xe = randn(size(F,1),vars); % E

% normalize the matrices
for f=1:3
    Xf{f} = sqrt(size(Xf{f},1))*Xf{f}/norm(Xf{f},'fro');
end
Xi = sqrt(size(Xi,1))*Xi/norm(Xi,'fro');
Xe = sqrt(size(Xe,1))*Xe/norm(Xe,'fro');

X = [];  
Xf1 = [];
Xf2 = [];
Xf3 = [];
Xi1 = [];
N = size(F,1); 
for i = 1:N     
    Xf1(i,:) = Xf{1}(F(i,1),:);
    Xf2(i,:) = Xf{2}(F(i,2),:);
    Xf3(i,:) = Xf{3}(F(i,1)+(F(i,3)-1)*length(levels{1}),:);
    Xi1(i,:) = Xi(F(i,1)+(F(i,2)-1)*length(levels{1}),:); 
end

X = k(1)*Xf1 + k(2)*Xf2 + k(3)*Xf3 + k(4)*Xi1 + k(5)*Xe; 

T = parglm_dist(X, F, {[1 2]},1,[],2,[],[],[],[1 3]);
legend('A','B','C(A)','AB')
axis([0.8 1.2 0 200])
saveas(gcf,'./Figures/Null3_1'); saveas(gcf,'./Figures/Null3_1.eps','epsc'); 

T.Source(2:5) = {'A','B','C(A)','AB'} 
table2latex(T,'./Figures/Table3_1.tex')

%% NULL distributions and ASCA table from data simulated with seed = 1. 
% The design matrix F contains a full factorial design with 
% four levels for A, three levels for B and four individuals in each cell 
% of C(A). Other inputs are M = 400, kA = kB = kC(A) = kAB = 0.1*theta 
% and kE = 1-theta, theta = 0, R = 1000, P = 200 and alpha = 0.05. The 
% number of levels of A, L_A, is duplicated.

clear
close all
clc

reps = 4;
vars = 400;
levels = {1:8,[1,2,3],1:reps};
theta = 0.5;
k = [theta*[.1,.1,.1,.1], (1-theta)]; 

F = create_design(levels,1);

rng(1);
    
for f=1:2 % main factors A and B
    Xf{f} = randn(length(levels{f}),vars);
end
Xf{3} = randn(length(levels{1})*length(levels{3}),vars);% C(A)
Xi = randn(length(levels{1})*length(levels{2}),vars); % AB
Xe = randn(size(F,1),vars); % E

% normalize the matrices
for f=1:3
    Xf{f} = sqrt(size(Xf{f},1))*Xf{f}/norm(Xf{f},'fro');
end
Xi = sqrt(size(Xi,1))*Xi/norm(Xi,'fro');
Xe = sqrt(size(Xe,1))*Xe/norm(Xe,'fro');

X = [];  
Xf1 = [];
Xf2 = [];
Xf3 = [];
Xi1 = [];
N = size(F,1); 
for i = 1:N     
    Xf1(i,:) = Xf{1}(F(i,1),:);
    Xf2(i,:) = Xf{2}(F(i,2),:);
    Xf3(i,:) = Xf{3}(F(i,1)+(F(i,3)-1)*length(levels{1}),:);
    Xi1(i,:) = Xi(F(i,1)+(F(i,2)-1)*length(levels{1}),:); 
end

X = k(1)*Xf1 + k(2)*Xf2 + k(3)*Xf3 + k(4)*Xi1 + k(5)*Xe; 

T = parglm_dist(X, F, {[1 2]},1,[],2,[],[],[],[1 3]);
legend('A','B','C(A)','AB')
axis([0.8 1.2 0 200])
saveas(gcf,'./Figures/Null3_2'); saveas(gcf,'./Figures/Null3_2.eps','epsc'); 

T.Source(2:5) = {'A','B','C(A)','AB'} 
table2latex(T,'./Figures/Table3_2.tex')

%% NULL distributions and ASCA table from data simulated with seed = 1. 
% The design matrix F contains a full factorial design with 
% four levels for A, three levels for B and four individuals in each cell 
% of C(A). Other inputs are M = 400, kA = kB = kC(A) = kAB = 0.1*theta 
% and kE = 1-theta, theta = 0, R = 1000, P = 200 and alpha = 0.05. The 
% number of levels of B, L_B, is duplicated.

clear
close all
clc

reps = 4;
vars = 400;
levels = {[1,2,3,4],1:6,1:reps};
theta = 0.5;
k = [theta*[.1,.1,.1,.1], (1-theta)]; 

F = create_design(levels,1);

rng(1);
    
for f=1:2 % main factors A and B
    Xf{f} = randn(length(levels{f}),vars);
end
Xf{3} = randn(length(levels{1})*length(levels{3}),vars);% C(A)
Xi = randn(length(levels{1})*length(levels{2}),vars); % AB
Xe = randn(size(F,1),vars); % E

% normalize the matrices
for f=1:3
    Xf{f} = sqrt(size(Xf{f},1))*Xf{f}/norm(Xf{f},'fro');
end
Xi = sqrt(size(Xi,1))*Xi/norm(Xi,'fro');
Xe = sqrt(size(Xe,1))*Xe/norm(Xe,'fro');

X = [];  
Xf1 = [];
Xf2 = [];
Xf3 = [];
Xi1 = [];
N = size(F,1); 
for i = 1:N     
    Xf1(i,:) = Xf{1}(F(i,1),:);
    Xf2(i,:) = Xf{2}(F(i,2),:);
    Xf3(i,:) = Xf{3}(F(i,1)+(F(i,3)-1)*length(levels{1}),:);
    Xi1(i,:) = Xi(F(i,1)+(F(i,2)-1)*length(levels{1}),:); 
end

X = k(1)*Xf1 + k(2)*Xf2 + k(3)*Xf3 + k(4)*Xi1 + k(5)*Xe; 

T = parglm_dist(X, F, {[1 2]},1,[],2,[],[],[],[1 3]);
legend('A','B','C(A)','AB')
axis([0.8 1.2 0 200])
saveas(gcf,'./Figures/Null3_3'); saveas(gcf,'./Figures/Null3_3.eps','epsc'); 

T.Source(2:5) = {'A','B','C(A)','AB'} 
table2latex(T,'./Figures/Table3_3.tex')

%% NULL distributions and ASCA table from data simulated with seed = 1. 
% The design matrix F contains a full factorial design with 
% four levels for A, three levels for B and four individuals in each cell 
% of C(A). Other inputs are M = 400, kA = kB = kC(A) = kAB = 0.1*theta 
% and kE = 1-theta, theta = 0, R = 1000, P = 200 and alpha = 0.05. The 
% number of replicates in C(A), r_{C(A)}, is triplicated.

clear
close all
clc

reps = 8;
vars = 400;
levels = {[1,2,3,4],[1,2,3],1:reps};
theta = 0.5;
k = [theta*[.1,.1,.1,.1], (1-theta)]; 

F = create_design(levels,1);

rng(1);
    
for f=1:2 % main factors A and B
    Xf{f} = randn(length(levels{f}),vars);
end
Xf{3} = randn(length(levels{1})*length(levels{3}),vars);% C(A)
Xi = randn(length(levels{1})*length(levels{2}),vars); % AB
Xe = randn(size(F,1),vars); % E

% normalize the matrices
for f=1:3
    Xf{f} = sqrt(size(Xf{f},1))*Xf{f}/norm(Xf{f},'fro');
end
Xi = sqrt(size(Xi,1))*Xi/norm(Xi,'fro');
Xe = sqrt(size(Xe,1))*Xe/norm(Xe,'fro');

X = [];  
Xf1 = [];
Xf2 = [];
Xf3 = [];
Xi1 = [];
N = size(F,1); 
for i = 1:N     
    Xf1(i,:) = Xf{1}(F(i,1),:);
    Xf2(i,:) = Xf{2}(F(i,2),:);
    Xf3(i,:) = Xf{3}(F(i,1)+(F(i,3)-1)*length(levels{1}),:);
    Xi1(i,:) = Xi(F(i,1)+(F(i,2)-1)*length(levels{1}),:); 
end

X = k(1)*Xf1 + k(2)*Xf2 + k(3)*Xf3 + k(4)*Xi1 + k(5)*Xe; 

T = parglm_dist(X, F, {[1 2]},1,[],2,[],[],[],[1 3]);
legend('A','B','C(A)','AB')
axis([0.8 1.2 0 200])
saveas(gcf,'./Figures/Null3_4'); saveas(gcf,'./Figures/Null3_4.eps','epsc'); 

T.Source(2:5) = {'A','B','C(A)','AB'} 
table2latex(T,'./Figures/Table3_4.tex')
