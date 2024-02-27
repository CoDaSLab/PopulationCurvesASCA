function [D, pD, powercurveo] = auxD(F, model, ordinal, coding, nested)

%
% INPUTS:
%
% F: [NxF] design matrix, cell or array, where columns correspond to 
% factors and rows to levels
%
% model: This paremeter is similar to 'model' of anovan. It could be:
%       'linear': only main effects are provided (by default)
%       'interaction': two order interactions are provided
%       'full': all potential interactions are provided
%       [1x1]: maximum order of interactions considered
%       [ix2]: array with two order interactions
%       cell: with each element a vector of factors
%
% ordinal: [1xF] whether factors are nominal or ordinal
%       0: nominal (default)
%       1: ordinal
% 
% coding: [1xF] type of coding of factors
%       0: sum/deviation coding (default)
%       1: reference coding (reference is the last level)
%
% nested: [nx2] pairs of neted factors, e.g., if factor 2 is nested in 1,
%   and 3 in 2, then nested = [1 2; 2 3]
%
%
% coded by: José Camacho (josecamacho@ugr.es)
% last modification: 20/Dec/23
%
% Copyright (C) 2023  Universidad de Granada
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

%% Arguments checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);

n_factors = size(F,2);                 % number of factors

if nargin < 2 || isempty(model), model = 'linear'; end;

if isequal(model,'linear')
    interactions = [];
end  
    
if isequal(model,'interaction')
    interactions = allinter(n_factors,2);
end    

if isequal(model,'full')
    interactions = allinter(n_factors,n_factors);
end    

if isnumeric(model) && isscalar(model) && model >= 2 && model <= n_factors
        interactions = allinter(n_factors,model);
end    

if isnumeric(model) && ~isscalar(model)
        interactions = {model};
end    

if iscell(model), interactions = model; end

if nargin < 3 || isempty(ordinal), ordinal = zeros(1,size(F,2)); end;
if nargin < 4 || isempty(coding), coding = zeros(1,size(F,2)); end;
if nargin < 5 || isempty(nested), nested = []; end;


%% Main code
                  
n_interactions      = length(interactions);      % number of interactions
n_factors           = size(F,2);                 % number of factors

% Create Design Matrix
n = 0;
D = [];
n = 1;
D = ones(size(F,1),1);
   

for f = 1 : n_factors
    if ordinal(f)
        D(:,n+1) = preprocess2D(F(:,f),1);
        powercurveo.factors{f}.Dvars = n+1;
        n = n + 1;
        powercurveo.factors{f}.order = 1;
    else
        if isempty(nested) || isempty(find(nested(:,2)==f)) % if not nested
            powercurveo.factors{f}.factors = [];
            uF = unique(F(:,f));
            powercurveo.n_levels(f) = length(uF);
            for i = 2:length(uF)
                D(find(ismember(F(:,f),uF(i))),n+i-1) = 1;
            end
            powercurveo.factors{f}.Dvars = n+(1:length(uF)-1);
            if coding(f) == 1
                D(find(ismember(F(:,f),uF(1))),powercurveo.factors{f}.Dvars) = 0;
            else
                D(find(ismember(F(:,f),uF(1))),powercurveo.factors{f}.Dvars) = -1;
            end
            n = n + length(uF) - 1;
            powercurveo.factors{f}.order = 1;
        else % if nested
            ind = find(nested(:,2)==f);
            ref = nested(ind,1);
            powercurveo.factors{f}.factors = [ref powercurveo.factors{ref}.factors];
            urF = unique(F(:,ref));
            powercurveo.n_levels(f) = 0;
            powercurveo.factors{f}.Dvars = [];
            for j = 1:length(urF)
                rind = find(ismember(F(:,ref),urF(j)));
                uF = unique(F(rind,f));
                powercurveo.n_levels(f) = powercurveo.n_levels(f) + length(uF);
                for i = 2:length(uF)
                    D(rind(find(ismember(F(rind,f),uF(i)))),n+i-1) = 1;
                end
                powercurveo.factors{f}.Dvars = [powercurveo.factors{f}.Dvars n+(1:length(uF)-1)];
                if coding(f) == 1
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = 0;
                else
                    D(rind(find(ismember(F(rind,f),uF(1)))),n+(1:length(uF)-1)) = -1;
                end
                n = n + length(uF) - 1;
            end   
            powercurveo.factors{f}.order = powercurveo.factors{ref}.order + 1;
        end
    end
end

for i = 1 : n_interactions
    Dout = computaDint(interactions{i},powercurveo.factors,D);
    D = [D Dout];
    powercurveo.interactions{i}.Dvars = n+1:size(D,2);
    powercurveo.interactions{i}.factors = interactions{i};
    n = size(D,2);
    powercurveo.interactions{i}.order = max(powercurveo.factors{interactions{i}(1)}.order,powercurveo.factors{interactions{i}(2)}.order) + 1;
    powercurveo.Dvars(powercurveo.interactions{i}.Dvars) = powercurveo.interactions{i}.order;
end

pD =  pinv(D'*D)*D';
    

end

%% Auxiliary function for interactions

function interactions = allinter(nF,order)
    
    if order > 2
        interactions = allinter(nF,order-1);
        for i = 1:length(interactions)
            for j = max(interactions{i})+1:nF
                interactions{end+1} = [interactions{i} j];
            end
        end
    else
        interactions = {};
        for i = 1:nF
            for j = i+1:nF
                interactions{end+1} = [i j];
            end
        end
    end
    
end
    
        
function Dout = computaDint(interactions,factors,D) % Compute coding matrix

    if length(interactions)>1
        deepD = computaDint(interactions(2:end),factors,D);
        Dout = [];
        for k = factors{interactions(1)}.Dvars
            for l = 1:size(deepD,2)
                Dout(:,end+1) = D(:,k).* deepD(:,l);
            end
        end
    else
        Dout = D(:,factors{interactions}.Dvars);
    end

end









