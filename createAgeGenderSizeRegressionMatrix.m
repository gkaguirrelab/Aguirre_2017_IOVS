function [ X, interactionX ] = createAgeGenderSizeRegressionMatrix(groupVector,polynomialDims,varargin)

numCovars=length(varargin);

% Create a balanced groupVector

positiveGroup=find(groupVector > 0);
numPositive=length(positiveGroup);
negativeGroup=find(groupVector < 0);
numNegative=length(negativeGroup);
groupVector(positiveGroup)=sum(groupVector(positiveGroup))/numPositive*(numPositive/numNegative);
groupVector(negativeGroup)=sum(groupVector(negativeGroup))/numNegative*(numNegative/numPositive);

% Create polynomial expansions of the base covariates
for c=1:numCovars
    baseCovar=cell2mat(varargin(c));
    subX(1,:)=baseCovar-mean(baseCovar);
    if polynomialDims(c)>1
        for p=2:polynomialDims(c)
            polyX=subX(1,:).^p;
            polyX=polyX-mean(polyX);
            polyX=polyX/max(polyX);
            [~,~,stats]=glmfit(subX',polyX');
            subX(p,:)=stats.resid;
        end
    end
    if c==1
        X=subX;
    else
        X=[X;subX];
    end
    clear subX
end

% Create a matrix to test for interactions
interactionX=X*0;

% Multiply main effect covariates by the group vector to
% create the interaction terms
for i=1:size(X,1)
    interactionX(i,:)=X(i,:).*groupVector';
end

% Render the interaction terms orthogonal to the main effect

for i=1:size(X,1)
    interaction=interactionX(i,:);
    mainEffect=X(i,:);
    [~,~,stats]=glmfit(mainEffect,interaction);
    interactionX(i,:)=stats.resid;
end

% detect any interaction terms that are all zeros and remove them

goodColumns=ones(size(X,1),1);
for i=1:size(X,1)
    if (max(interactionX(i,:))-min(interactionX(i,:)))<0.01
        goodColumns(i)=0;
    end
end

interactionX=interactionX(find(goodColumns==1),:);

gribble=1;
end

