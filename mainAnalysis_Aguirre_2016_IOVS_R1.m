% Figures and Analyses for LCA1  paper

clear;
close all;

%  CAN SWITCH BETWEEN PROCESSING THE AUDITORY AND VISUAL DATA AT THIS STAGE

figureNameStem='LCA1_Fig';

data_type = 'visual';
%data_type = 'auditory';

% Place to save figures, ideally in the DropBox hierarchy
OutFileStem=[''];
FigIndex=1;

% The test / retest data are for the visual data set only

indextest=linspace(1,24,24);
indexretest=linspace(25,48,24);

% Create a header to mark out the text output

TAB=char(9);  % Set up a tab to be used in text output
fprintf('\n\n\n***************************************************\n\n')

% Load a .csv data file that has all the integrated subject data and
% information. The measures in this file are the product of scripts
% produced by Rito Datta that act upon the anatomical and fMRI data.

PrimaryDataFileName=['Blindness_DataMatrix.csv'];
T = readtable(PrimaryDataFileName);
clear filename;

% Extract a vector that indicates if the subject was studied under
% Protocol:
% 1 - prior to June 23, 2008, or for the "weather" or "migraine" subjects
% 2 - July 2008 or later
% 3 - MPRAGE file from Oxford

protocolNumber = T.ProtocolNumber;

% Extract a set of vectors from this table to idenitify subjects by group,
% gender, age, etc.

gender = double(strcmp(T.Gender,'Male')); % Female = 0, Male = 1
ages = T.Age;
blindness_onset_max_years=[T.VisionLossOnset T.AgeAtMaxLoss];
blindnessOnset=T.VisionLossOnset;
SubjectID=T.ExternalSubjectID;
NumSubjects = length(SubjectID);

lca1PlotMarkers={'^','d','v','h','s','o'};

% This set of variables contain indicies that identify sub-populations of
% the total set of subjects studied

indexsight=find(strcmp(T.Group,'Sighted')==1)';
indexanophthalmic=find(strcmp(T.Group,'Anophthalmic')==1)';
indexcongenital=find(strcmp(T.Group,'Congenital')==1)';
indexpostnatal=find(strcmp(T.Group,'Postnatal')==1)';

indexrpe65=find(strcmp(T.Group,'rpe65')==1)';
indexcrb1=find(strcmp(T.Group,'crb1')==1)';
indexcep290=find(strcmp(T.Group,'CEP290')==1)';
indexlca1=find(strcmp(T.Group,'LCA1')==1)'; % note that the LCA1 subjects are not included in the list of general LCA subjects

indexlca=[indexrpe65 indexcrb1 indexcep290];
indexblind=[indexanophthalmic indexcongenital indexpostnatal indexlca];
indexblindWithLCA1=[indexanophthalmic indexcongenital indexpostnatal indexlca indexlca1];

indexEarlyBlind=[indexanophthalmic indexcongenital indexlca];

indexpostnatal_beforepuberty = find( (strcmp(T.Group,'Postnatal').*(T.VisionLossOnset<=14)) ==1 )';
indexpostnatal_afterpuberty = find( (strcmp(T.Group,'Postnatal').*(T.VisionLossOnset>=15)) ==1 )';

% Identify the sighted subjects studied under protocol 1 and 2
indexsightProtocol1=find((strcmp(T.Group,'Sighted')==1).*(protocolNumber==1))';
indexsightProtocol2=find((strcmp(T.Group,'Sighted')==1).*(protocolNumber==2))';



% This is a matrix of the non-MPRAGE derived data

othermeasures=NaN([NumSubjects,3]);

othermeasures(:,1)=T.V1meanCBF;       % V1 CBF values already scaled by global CBF
othermeasures(:,2)=(T.lhV1SenNoiseBeta + T.rhV1SenNoiseBeta + T.lhV1RevNoiseBeta + T.rhV1RevNoiseBeta);
othermeasures(:,3)=(T.FASplenium+T.FAOpticRadiation)/2;
%othermeasures(:,3)=(T.FAOpticRadiation);
%othermeasures(:,3)=(T.FASplenium);
othermeasures(:,4)=T.restingHierarchical-T.restingHomotopic;
othermeasures(:,5)=T.restingHierarchical;
othermeasures(:,6)=T.restingHomotopic;

% Now assemble a matrix that contains the core anatomical values,
%  prior to adjustments for age, gender, overall size

if strcmp(data_type,'visual')
    raw_data=zeros(NumSubjects,8);
    measures={'lhV1thick','rhV1thick','lhV1SurfArea','rhV1SurfArea','LHpericalVol','RHpericalVol','ChiasmVol','LGNjacobian'};%,'SpleniumVol'};
    measure_type={'thickness','thickness','area','area','volume','volume','volume','jacobian'};%,'volume'};
    raw_data(:,1)= T.lhV1thick;
    raw_data(:,2)= T.rhV1thick;
    raw_data(:,3)= T.lhV1SurfArea;
    raw_data(:,4)= T.rhV1SurfArea;
    raw_data(:,5)= T.LHpericalVol;
    raw_data(:,6)= T.RHpericalVol;
    raw_data(:,7)= T.ChiasmVol;
    raw_data(:,8)= T.LGNjacobian;
    %raw_data(:,9)= T.SpleniumVol;
    
    fprintf('\n\nVISUAL DATA ANALYSIS\n\n');
    OutFileStem=[OutFileStem 'Figures/'];
end

if strcmp(data_type,'auditory')
    raw_data=zeros(NumSubjects,7);
    measures={'lhA1thick','rhA1thick','lhA1SurfArea','rhA1SurfArea','LHtranstemporalVol','RHtranstemporalVol','MGNjacobian'};
    measure_type={'thickness' 'thickness' 'area' 'area' 'volume' 'volume' 'jacobian'};
    raw_data(:,1)= T.lhA1thick;
    raw_data(:,2)= T.rhA1thick;
    raw_data(:,3)= T.lhA1SurfArea;
    raw_data(:,4)= T.rhA1SurfArea;
    raw_data(:,5)= T.LHtranstemporalVol;
    raw_data(:,6)= T.RHtranstemporalVol;
    raw_data(:,7)= T.MGNjacobian;
    
    fprintf('\n\nAUDITORY DATA ANALYSIS\n\n');
    OutFileStem=[OutFileStem 'Figures/'];
end

sizer=size(raw_data);
NumMeasures=sizer(2);
clear sizer;


% Obtain some whole brain structural measures that will be used to adjust
% the data for individual differences in size

scalerNames={'global thickness','global area','supratentorial volume','intracranial volume'};
thicknessScaler=(T.lhCortexthick+T.rhCortexthick)/2;  % Scale by whole brain average cortical thickness
areaScaler=T.lhCortexSurfArea+T.rhCortexSurfArea; % Scale by whole brian surface area
volumeScalerA=T.SupraTenVol;
volumeScalerB=T.EstimatedTotalIntraCranialVol;

clear T;  % Done with the table


%% Create a little subject demographics table

fprintf('\n\nSubject Demographics:\n\n');
fprintf(['Sighted: ' num2str(sum(gender(indexsight))) 'm / ' num2str(length(indexsight)-sum(gender(indexsight))) 'f, age= ' num2str(mean(ages(indexsight)),'%1.0f') ' ± ' num2str(std(ages(indexsight)),'%1.0f') '\n']);
fprintf(['LCA1: ' num2str(sum(gender(indexlca1))) 'm / ' num2str(length(indexlca1)-sum(gender(indexlca1))) 'f, age= ' num2str(mean(ages(indexlca1)),'%1.0f') ' ± ' num2str(std(ages(indexlca1)),'%1.0f') '\n']);
fprintf(['All blind (not LCA1): ' num2str(sum(gender(indexblind))) 'm / ' num2str(length(indexblind)-sum(gender(indexblind))) 'f, age= ' num2str(mean(ages(indexblind)),'%1.0f') ' ± ' num2str(std(ages(indexblind)),'%1.0f') '\n']);
fprintf(['  [Anophthalmic, congenital, LCA]: ' num2str(sum(gender(indexEarlyBlind))) 'm / ' num2str(length(indexEarlyBlind)-sum(gender(indexEarlyBlind))) 'f, age= ' num2str(mean(ages(indexEarlyBlind)),'%1.0f') ' ± ' num2str(std(ages(indexEarlyBlind)),'%1.0f') '\n']);
fprintf(['  Postnatal: ' num2str(sum(gender(indexpostnatal))) 'm / ' num2str(length(indexpostnatal)-sum(gender(indexpostnatal))) 'f, age= ' num2str(mean(ages(indexpostnatal)),'%1.0f') ' ± ' num2str(std(ages(indexpostnatal)),'%1.0f') '\n']);
fprintf('\n\n');



%% Adjust the data to account for age, gender, and whole brain scaling effects
% Store and report the beta values associated with these adjustments.

% Create a vector that models the difference between protocol 1 and 2 at
% the Penn site

protocolCovariate=zeros(NumSubjects,1);
protocolCovariate(indexsightProtocol1)=1/length(indexsightProtocol1);
protocolCovariate(indexsightProtocol2)=(-1)/length(indexsightProtocol2);


adjusted_data=zeros(NumSubjects,NumMeasures);

indexAdjustGroup=[indexblindWithLCA1 indexsight];
groupVector=zeros(NumSubjects,1);
groupVector(indexblindWithLCA1)=-1;
groupVector(indexsight)=1;
%fprintf('_Effects of age and gender in the scaled BLIND and SIGHTED data_\n');
%fprintf('Note that the effect of protcol was only modeled in the sighted control subjects\n');
%fprintf('Beta (p), for each covariate [age, age^2, age^3, protocol1or2, gender, sizescalerA, sizescalerA^2, sizescalerA^3, sizerscalerB, sizescalerB^2, sizerscalerB^3]\n\n');

for i=1:NumMeasures
    if strcmp(measure_type(i),'thickness')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup), gender(indexAdjustGroup), thicknessScaler(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'area')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), areaScaler(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'volume')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), volumeScalerA(indexAdjustGroup), volumeScalerB(indexAdjustGroup));
    end
    if strcmp(measure_type(i),'jacobian')
        [regressionMatrix,~]=createAgeGenderSizeRegressionMatrix(groupVector,[3,1,1,3,3],ages(indexAdjustGroup),protocolCovariate(indexAdjustGroup),gender(indexAdjustGroup), volumeScalerA(indexAdjustGroup), volumeScalerB(indexAdjustGroup));
    end
    regressionMatrix=[regressionMatrix;ones(1,length(indexAdjustGroup))];
    y=raw_data(indexAdjustGroup,i);
    weights=y*0+1;
    weights(indexlca1)=0; % set the weights on the LCA1 subjects to zero, so they do not influence the calculation of beta values.
    [b,~,stats] = glmfit(regressionMatrix',y,'normal','constant','off','weights',weights);
    adjusted_data((indexAdjustGroup),i)=stats.resid+b(end);
    Outline=[char(measures(i)) TAB];
    for betas=1:length(b)-1
        Outline=[Outline num2str(b(betas),'%1.2e') ' (' num2str(stats.p(betas),'%1.3f') ')' TAB];
    end
%    fprintf([Outline  '\n']);
    clear b;
    clear stats;
    clear y;
    clear regressionMatrix;
end
fprintf('\n\n');



clear raw_data;  % done with the raw measures
clear indexgroup;
clear groups;
clear scaled_data;

data=adjusted_data;
clear adjusted_data;

%% Average the data across hemispheres

numClust=5;
clustNames={'V1thick','V1SurfArea','pericalVol','ChiasmVol','LGNjacobian'};

all_score=NaN([NumSubjects numClust]);
all_score(:,1)=mean(data(:,[1 2]),2);
all_score(:,2)=mean(data(:,[3 4]),2);
all_score(:,3)=mean(data(:,[5 6]),2);
all_score(:,4)=mean(data(:,[7]),2);
all_score(:,5)=mean(data(:,[8]),2);



% Report Kruskall-Wallis tests

group=all_score(:,1)*nan;
group(indexsight)=1;
group(indexEarlyBlind)=2;
group(indexlca1)=3;

fprintf('Kruskall-Wallis test for an overall effect of group for each cluster:\n');
for i=1:numClust
    [p,tbl,stats] = kruskalwallis(all_score(:,i),group,'off');
    fprintf(['measure ' clustNames{i} ': Chi-squared (2) = ' num2str(tbl{2,5},'%1.5f') ', p = ' num2str(tbl{2,6},'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('follow-up test for blind different from sighted:\n');
for i=1:numClust
    [p,tbl,stats] = kruskalwallis(all_score(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' clustNames{i} ': p= ' num2str(c(1,6),'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('follow-up test for lca1 different from sighted:\n');
for i=1:numClust
    [p,tbl,stats] = kruskalwallis(all_score(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' clustNames{i} ': p= ' num2str(c(2,6),'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('follow-up test for lca1 different from blind:\n');
for i=1:numClust
    [p,tbl,stats] = kruskalwallis(all_score(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' clustNames{i} ': p= ' num2str(c(3,6),'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');


% Create bar plots of the group means and SEM for each cluster

figtmp=figure('name','Group Means ± SEM each dimension');
titles=clustNames;

PCLabels={'sighted' 'blind-all' 'lca1'};

group=all_score(:,1)*nan;
group(indexsight)=1;
group(indexEarlyBlind)=2;
group(indexlca1)=3;

ClusterPlotLocation=[6,5,4,2,3];

yLims = [ [ 0 500 ]; [0 1.5]; [0 5000]; [0 3000]; [0 2.5] ];

for i=1:numClust
    subplot(2,3,ClusterPlotLocation(i));
    % calculate a jitter for the plot points, based upon their density
    [histN,histEdges]=histcounts(all_score(indexsight,i),10);
    j=all_score(indexsight,i)*0;
    for b=1:10
        idx=find( (all_score(indexsight,i)>=histEdges(b)) .* (all_score(indexsight,i)<histEdges(b+1)));
        if ~isempty(idx)
            maxShift=(histN(b)/max(histN)) * 0.15;
            jitters= -maxShift : (maxShift*2)/length(idx) : maxShift;
            [~,sortIdx]=sort(abs(jitters));
            jitters=jitters(sortIdx);
            jitters=jitters(1:length(idx));
            j(idx)= jitters;
        end
    end
    plot(.75+j,all_score(indexsight,i),'Marker','.','LineStyle','none','MarkerEdgeColor',[.75 .75 .75],'MarkerFaceColor','none')
    hold on
    
        % calculate a jitter for the plot points, based upon their density
    [histN,histEdges]=histcounts(all_score(indexEarlyBlind,i),10);
    j=all_score(indexEarlyBlind,i)*0;
    for b=1:10
        idx=find( (all_score(indexEarlyBlind,i)>=histEdges(b)) .* (all_score(indexEarlyBlind,i)<histEdges(b+1)));
        if ~isempty(idx)
            maxShift=(histN(b)/max(histN)) * 0.1;
            jitters= -maxShift : (maxShift*2)/length(idx) : maxShift;
            [~,sortIdx]=sort(abs(jitters));
            jitters=jitters(sortIdx);
            jitters=jitters(1:length(idx));
            j(idx)= jitters;
        end
    end

    plot(1.75+j,all_score(indexEarlyBlind,i),'Marker','.','LineStyle','none','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor','none')
    % Add plot symbols for each LCA1 subject
    for n=1:6
       subName= SubjectID{indexlca1(n)};
       markerIndex=str2num(subName(end));
       plot(2.75,all_score(indexlca1(n),i),'Marker',lca1PlotMarkers{markerIndex},'MarkerEdgeColor',[.75 0 0],'MarkerFaceColor','none')
    end

    boxplot(all_score(:,i),group,'Notch','on','MedianStyle','line','BoxStyle','filled','Symbol','','Whisker',0,'Colors',[[.5 .5 .5];[.25 .25 .25];[1 0 0]])
ylim(yLims(ClusterPlotLocation(i)-1,:));
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([1 1.5 1])
    box off;
    title([titles{i}]);
    hold off
    
end

saveas(figtmp, [OutFileStem figureNameStem num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


%%%%%%%%%%%%%%%%%%%%%%%
%% Time to work on the "other measures"
%%%%%%%%%%%%%%%%%%%%%%%


%% Report the number of subjects who have each kind of "other" measure

fprintf('\n\nNumber of subjects with additional brain measures:\n\n');
fprintf( '         CBF  Sentences FA Rest\n');
fprintf(['Sighted: ' num2str(length(indexsight)-sum(isnan(othermeasures(indexsight,:)))) '\n']);
fprintf(['Blind-all:   ' num2str(length(indexblind)-sum(isnan(othermeasures(indexblind,:)))) '\n']);
fprintf(['Blind-early:   ' num2str(length(indexEarlyBlind)-sum(isnan(othermeasures(indexEarlyBlind,:)))) '\n']);
fprintf(['LCA1:   ' num2str(length(indexlca1)-sum(isnan(othermeasures(indexlca1,:)))) '\n']);
fprintf('\n\n');


%% Report the means, and t-test of the sighted vs. the blind

% Obtain the mean and SEM of the non-MPRAGE measures, broken down
%  by sub-group. Create a plot of this.

LabelsMeasures={'CBF' 'Sentences' 'FA' 'Rest_Hier-Homo' 'Rest_Hier' 'Rest_Homo'};
LabelsGroups={'Sight' 'Blind-all' 'LCA1'};

% Use a common ranges to plot these results
 
YRanges=[0 3; -5 10; 0.3 0.6];
XRanges=YRanges;


figtmp=figure('name','Median of other measures by group ±IQR');
titles={'CBF' 'Sentences' 'FA' 'Rest_Hier-Homo' 'Rest_Hier' 'Rest_Homo'};
yLims = [ [ 0 2.5 ]; [-2.5 10]; [0 0.6]; [-1 1]; [0 1.5]; [0 1.5] ];

for i=1:6
    subplot(2,3,i);
    % calculate a jitter for the plot points, based upon their density
    [histN,histEdges]=histcounts(othermeasures(indexsight,i),10);
    j=othermeasures(indexsight,i)*0;
    for b=1:10
        idx=find( (othermeasures(indexsight,i)>=histEdges(b)) .* (othermeasures(indexsight,i)<histEdges(b+1)));
        if ~isempty(idx)
            maxShift=(histN(b)/max(histN)) * 0.15;
            jitters= -maxShift : (maxShift*2)/length(idx) : maxShift;
            [~,sortIdx]=sort(abs(jitters));
            jitters=jitters(sortIdx);
            jitters=jitters(1:length(idx));
            j(idx)= jitters;
        end
    end
    plot(.75+j,othermeasures(indexsight,i),'Marker','.','LineStyle','none','MarkerEdgeColor',[.75 .75 .75],'MarkerFaceColor','none')
    hold on
    
        % calculate a jitter for the plot points, based upon their density
    [histN,histEdges]=histcounts(othermeasures(indexEarlyBlind,i),10);
    j=othermeasures(indexEarlyBlind,i)*0;
    for b=1:10
        idx=find( (othermeasures(indexEarlyBlind,i)>=histEdges(b)) .* (othermeasures(indexEarlyBlind,i)<histEdges(b+1)));
        if ~isempty(idx)
            maxShift=(histN(b)/max(histN)) * 0.1;
            jitters= -maxShift : (maxShift*2)/length(idx) : maxShift;
            [~,sortIdx]=sort(abs(jitters));
            jitters=jitters(sortIdx);
            jitters=jitters(1:length(idx));
            j(idx)= jitters;
        end
    end

    plot(1.75+j,othermeasures(indexEarlyBlind,i),'Marker','.','LineStyle','none','MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor','none')
    % Add plot symbols for each LCA1 subject
    for n=1:6
       subName= SubjectID{indexlca1(n)};
       markerIndex=str2num(subName(end));
       plot(2.75,othermeasures(indexlca1(n),i),'Marker',lca1PlotMarkers{markerIndex},'MarkerEdgeColor',[.75 0 0],'MarkerFaceColor','none')
    end

    boxplot(othermeasures(:,i),group,'Notch','on','MedianStyle','line','BoxStyle','filled','Symbol','','Whisker',0,'Colors',[[.5 .5 .5];[0 0 0];[1 0 0]])
ylim(yLims(i,:));
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([1 1.5 1])
    box off;
    title([titles{i}]);
    hold off
end

saveas(figtmp, [OutFileStem figureNameStem num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Report Kruskall-Wallis tests

group=othermeasures(:,1)*nan;
group(indexsight)=1;
group(indexEarlyBlind)=2;
group(indexlca1)=3;

fprintf('Kruskall-Wallis test for an overall effect of group for each cluster:\n');
for i=1:4
    [p,tbl,stats] = kruskalwallis(othermeasures(:,i),group,'off');
    fprintf(['measure ' LabelsMeasures{i} ': Chi-squared (2) = ' num2str(tbl{2,5},'%1.5f') ', p = ' num2str(tbl{2,6},'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('follow-up test for blind different from sighted:\n');
for i=1:4
    [p,tbl,stats] = kruskalwallis(othermeasures(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(1,6),'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('follow-up test for lca1 different from sighted:\n');
for i=1:4
    [p,tbl,stats] = kruskalwallis(othermeasures(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(2,6),'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('follow-up test for lca1 different from blind:\n');
for i=1:4
    [p,tbl,stats] = kruskalwallis(othermeasures(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(3,6),'%1.7f') '\n']);
    clear d;
end
fprintf('\n\n');





% Check that the FA value comparisons limited to just those subjects
% studied with protocol 2

fprintf('p value of Kruskall-Wallis tests for just those subjects with protocol 2 data:\n');
for i=3:3
    
    
tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexEarlyBlind)=1;
indexEarlyBlind_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexsight)=1;
indexsight_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexlca1)=1;
indexlca1_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);


group=othermeasures(:,1)*nan;
group(indexsight_protocol2)=1;
group(indexEarlyBlind_protocol2)=2;
group(indexlca1_protocol2)=3;

fprintf('follow-up test for lca1 different from sighted:\n');

    [p,tbl,stats] = kruskalwallis(othermeasures(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(2,6),'%1.7f') '\n']);
    clear d;

fprintf('\n\n');

fprintf('follow-up test for lca1 different from blind:\n');

    [p,tbl,stats] = kruskalwallis(othermeasures(:,i),group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(3,6),'%1.7f') '\n']);
    clear d;

fprintf('\n\n');

end
fprintf('\n\n');


%% Conduct an analysis in which we account for effects of age in the FA
% data, and test if the same group effects hold.

fprintf('Regress out the effect of age upon FA values and repeat the test:\n');

faData=othermeasures(:,3);
agevector=ages-sum(ages)/length(ages);
[~,~,faDataAgeRegressed,~] = regress(faData,agevector);


% Report t-tests comparing the groups in the age-regressed FA values

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexEarlyBlind)=1;
indexEarlyBlind_protocol2=find(~isnan(faDataAgeRegressed).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexsight)=1;
indexsight_protocol2=find(~isnan(faDataAgeRegressed).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexlca1)=1;
indexlca1_protocol2=find(~isnan(faDataAgeRegressed).*(protocolNumber==2).*tmpSubjectArray);

group=othermeasures(:,1)*nan;
group(indexsight_protocol2)=1;
group(indexEarlyBlind_protocol2)=2;
group(indexlca1_protocol2)=3;

fprintf('follow-up test for lca1 different from sighted:\n');

    [p,tbl,stats] = kruskalwallis(faDataAgeRegressed,group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(2,6),'%1.7f') '\n']);
    clear d;

fprintf('\n\n');

fprintf('follow-up test for lca1 different from blind:\n');

    [p,tbl,stats] = kruskalwallis(faDataAgeRegressed,group,'off');
    c = multcompare(stats,'Display','off');
    fprintf(['measure ' LabelsMeasures{i} ': p= ' num2str(c(3,6),'%1.7f') '\n']);
    clear d;

fprintf('\n\n');