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
fprintf(['  [Anophthalmic, congenital, LCA]: ' num2str(sum(gender(indexblind))) 'm / ' num2str(length(indexblind)-sum(gender(indexblind))) 'f, age= ' num2str(mean(ages(indexblind)),'%1.0f') ' ± ' num2str(std(ages(indexblind)),'%1.0f') '\n']);
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




%% Create bar plots of the group means and SEM for measure

figtmp=figure('name','Group Means ± SEM each measure');

groupLabels={'sighted' 'blind-all' 'LCA1'};


for i=1:NumMeasures
    subplot(3,3,i);
    tmpGroupDimScoreMeans=[mean(adjusted_data(indexsight,i)) mean(adjusted_data(indexblind,i)) mean(adjusted_data(indexlca1,i))];
    tmpGroupDimScoreSEM=[std(adjusted_data(indexsight,i))/sqrt(length(indexsight)) std(adjusted_data(indexblind,i))/sqrt(length(indexblind)) std(adjusted_data(indexlca1,i))/sqrt(length(indexlca1))];
    bar([1 2 3],tmpGroupDimScoreMeans);
    hold on
    errorbar([1 2 3],tmpGroupDimScoreMeans,tmpGroupDimScoreSEM,'o');
    hold on
    xlim([0 4]);
    set(gca,'Xtick',1:3,'XTickLabel',groupLabels);
    pbaspect([2 1 1])
    box off;
    title([measures{i}]);
end

saveas(figtmp, [OutFileStem figureNameStem num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


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


%% Report the mean and SEM of the measures in the sighted, blind-all, blind-early, and LCA1 groups,

fprintf('Mean ± SEM for each adjusted measure for each group:\n');
fprintf('        Sighted     blind-all    blind-early    LCA1:\n');
for i=1:numClust
    Outline=['Measure ' clustNames{i} '  ' ];
    Outline=[Outline num2str(mean(all_score(indexsight,i)),'%1.2f') '±' num2str(std(all_score(indexsight,i))/sqrt(length(indexsight)),'%1.2f') '  '];
    Outline=[Outline num2str(mean(all_score(indexblind,i)),'%1.2f') '±' num2str(std(all_score(indexblind,i))/sqrt(length(indexblind)),'%1.2f') '  '];
    Outline=[Outline num2str(mean(all_score(indexEarlyBlind,i)),'%1.2f') '±' num2str(std(all_score(indexEarlyBlind,i))/sqrt(length(indexEarlyBlind)),'%1.2f') '  '];
    Outline=[Outline num2str(mean(all_score(indexlca1,i)),'%1.2f') '±' num2str(std(all_score(indexlca1,i))/sqrt(length(indexlca1)),'%1.2f') '  '];
    fprintf([Outline '\n']);
end



% Report t-tests comparing the groups in the clustered values

fprintf('p value of a t-test between the sighted and LCA1 for each cluster:\n');
for i=1:numClust
    [~,pVal]=ttest2(all_score(indexsight,i),all_score(indexlca1,i));
    fprintf(['measure ' clustNames{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('p value of a t-test between all the blind and LCA1 for each cluster:\n');
for i=1:numClust
    [~,pVal]=ttest2(all_score(indexblind,i),all_score(indexlca1,i));
    fprintf(['measure ' clustNames{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
fprintf('\n\n');

fprintf('p value of a t-test between the congenitally blind and LCA1 for each cluster:\n');
for i=1:numClust
    [~,pVal]=ttest2(all_score(indexEarlyBlind,i),all_score(indexlca1,i));
    fprintf(['measure ' clustNames{i} ': ' num2str(pVal,'%1.5f') '\n']);
    clear d;
end
fprintf('\n\n');


% Create bar plots of the group means and SEM for each cluster

figtmp=figure('name','Group Means ± SEM each dimension');
titles=clustNames;

PCLabels={'sighted' 'blind-all' 'lca1'};

for i=1:numClust
    subplot(3,2,i);
    tmpGroupDimScoreMeans=[mean(all_score(indexsight,i)) mean(all_score(indexblind,i)) mean(all_score(indexlca1,i))];
    tmpGroupDimScoreSTD=[std(all_score(indexsight,i)) std(all_score(indexblind,i)) std(all_score(indexlca1,i)) ];
    tmpGroupDimScoreSEM=[std(all_score(indexsight,i))/sqrt(length(indexsight)) std(all_score(indexblind,i))/sqrt(length(indexblind)) std(all_score(indexlca1,i))/sqrt(length(indexlca1)) ];
    bar([1 2 3],tmpGroupDimScoreMeans,'EdgeColor','none','FaceColor',[0.5 .5 .5]);
    hold on
    errorbar([1 2 3],tmpGroupDimScoreMeans,tmpGroupDimScoreSEM,'o');
    hold on
    xlim([0 6]);
    set(gca,'Xtick',1:3,'XTickLabel',PCLabels);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
    
    % Add plot symbols for each LCA1 subject
    for n=1:6
       subName= SubjectID{indexlca1(n)};
       markerIndex=str2num(subName(end));
       plot(4,all_score(indexlca1(n),i),'Marker',lca1PlotMarkers{markerIndex},'MarkerEdgeColor','red','MarkerFaceColor','none')
    end
    
end

saveas(figtmp, [OutFileStem figureNameStem num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


%% Test if the mean of the z-trasformed measures of the LGN, pericalcarine
% volume, and V1 surface area differ between the LCA1 and the blind

% Z-transform the data within each measure
%  The adjustment is performed on all subjects (including LCA1), but the
%  parameters of the adjustment are calculated excluding the LCA1 subjects.

z_all_score=NaN([NumSubjects numClust]);

for i=1:numClust
    z_all_score(:,i)=(all_score(:,i) - mean(all_score([indexblind indexsight],i))) / std(all_score([indexblind indexsight],i));
end

% Report a t-test 

fprintf('p value of a t-test between all blind and LCA1 for average of z-scored LGN, pericalc, and V1 area:\n');
    [~,p,~,stats]=ttest2(mean(z_all_score(indexlca1,2:4),2),mean(z_all_score(indexblind,2:4),2));
    fprintf(['t-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.3f') '\n']);
fprintf('\n\n');

fprintf('p value of a t-test between early blind and LCA1 for average of z-scored LGN, pericalc, and V1 area:\n');
    [~,p,~,stats]=ttest2(mean(z_all_score(indexlca1,2:4),2),mean(z_all_score(indexEarlyBlind,2:4),2));
    fprintf(['t-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.3f') '\n']);
fprintf('\n\n');



%%%%%%%%%%%%%%%%%%%%%%%
%% Time to work on the "other measures"
%%%%%%%%%%%%%%%%%%%%%%%


%% Report the number of subjects who have each kind of "other" measure

fprintf('\n\nNumber of subjects with additional brain measures:\n\n');
fprintf( '         CBF  Sentences FA\n');
fprintf(['Sighted: ' num2str(length(indexsight)-sum(isnan(othermeasures(indexsight,:)))) '\n']);
fprintf(['Blind-all:   ' num2str(length(indexblind)-sum(isnan(othermeasures(indexblind,:)))) '\n']);
fprintf(['Blind-early:   ' num2str(length(indexblind)-sum(isnan(othermeasures(indexEarlyBlind,:)))) '\n']);
fprintf(['LCA1:   ' num2str(length(indexlca1)-sum(isnan(othermeasures(indexlca1,:)))) '\n']);
fprintf('\n\n');


%% Report the means, and t-test of the sighted vs. the blind

% Obtain the mean and SEM of the non-MPRAGE measures, broken down
%  by sub-group. Create a plot of this.

LabelsMeasures={'CBF' 'Sentences' 'FA'};
LabelsGroups={'Sight' 'Blind-all' 'LCA1'};

OtherMeasureMean=zeros(3,3);
OtherMeasureSEM=zeros(3,3);

% Use a common ranges to plot these results
 
YRanges=[0 3; -5 10; 0.3 0.6];
XRanges=YRanges;

for i=1:3 % Loop across measures
    OtherMeasureMean(i,1)=nanmean(othermeasures(indexsight,i));
    OtherMeasureSEM(i,1)=nanstd(othermeasures(indexsight,i)) / sqrt(length(indexsight));
    
    OtherMeasureMean(i,2)=nanmean(othermeasures(indexblind,i));
    OtherMeasureSEM(i,2)=nanstd(othermeasures(indexblind,i)) / sqrt(length(indexblind));
    
    OtherMeasureMean(i,3)=nanmean(othermeasures(indexlca1,i));
    OtherMeasureSEM(i,3)=nanstd(othermeasures(indexlca1,i)) / sqrt(length(indexlca1));
end

figtmp=figure('name','Mean of other measures by group ±SEM');
titles={'CBF' 'Sentences' 'FA'};

for i=1:3
    subplot(3,1,i);
    
    bar([1 2 3],OtherMeasureMean(i,:),'EdgeColor','none','FaceColor',[0.5 .5 .5]);
    hold on
    errorbar([1 2 3],OtherMeasureMean(i,:),OtherMeasureSEM(i,:),'o');
    hold on
    xlim([0 5]);
    set(gca,'Xtick',1:3,'XTickLabel',LabelsGroups);
    pbaspect([2 1 1])
    box off;
    title([titles{i}]);
    
        % Add plot symbols for each LCA1 subject
    for n=1:6
       subName= SubjectID{indexlca1(n)};
       markerIndex=str2num(subName(end));
       plot(4,othermeasures(indexlca1(n),i),'Marker',lca1PlotMarkers{markerIndex},'MarkerEdgeColor','red','MarkerFaceColor','none')
    end
    
end

saveas(figtmp, [OutFileStem figureNameStem num2str(FigIndex) '.pdf'], 'pdf');
FigIndex=FigIndex+1;


% Report t-tests comparing the groups in the other measure

fprintf('p value of a t-test between the all-blind and sighted for each other measure:\n');
for i=1:3
    [~,p,~,stats]=ttest2(othermeasures(indexblind,i),othermeasures(indexsight,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.4f') '\n']);
end
fprintf('\n\n');


fprintf('p value of a t-test between the sighted and LCA1 for each other measure:\n');
for i=1:3
    [~,p,~,stats]=ttest2(othermeasures(indexsight,i),othermeasures(indexlca1,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.4f') '\n']);
end
fprintf('\n\n');

fprintf('p value of a t-test between the all-blind and LCA1 for each other measure:\n');
for i=1:3
    [~,p,~,stats]=ttest2(othermeasures(indexblind,i),othermeasures(indexlca1,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.4f') '\n']);
end
fprintf('\n\n');

fprintf('p value of a t-test between the early-blind and LCA1 for each other measure:\n');
for i=1:3
    [~,p,~,stats]=ttest2(othermeasures(indexEarlyBlind,i),othermeasures(indexlca1,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.4f') '\n']);
end
fprintf('\n\n');


% Make a little summary data table for the LCA1 subjects
fprintf('Data table for the GUCY2D-LCA subjects:\n');
fprintf('       V1 thickness   Optic chiasm volume    Cross Modal response    FA \n');
for i=1:length(indexlca1)
    Outline=SubjectID{indexlca1(i)};
    fprintf([Outline ': ' num2str(all_score(indexlca1(i),1),'%1.2f') ' & ' num2str(all_score(indexlca1(i),4),'%1.0f') ' & ' num2str(othermeasures(indexlca1(i),2),'%1.2f') ' & ' num2str(othermeasures(indexlca1(i),3),'%1.2f') '\n']);
end
fprintf('\n\n');

% Check that the FA value comparisons limited to just those subjects
% studied with protocol 2

fprintf('p value of a t-test between the blind and sighted for each other measure for just those subjects with protocol 2 data:\n');
for i=1:3
    
    
tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexblind)=1;
indexblind_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexsight)=1;
indexsight_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexlca1)=1;
indexlca1_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);


    [~,p,~,stats]=ttest2(othermeasures(indexblind_protocol2,i),othermeasures(indexsight_protocol2,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.3f') '\n']);
end
fprintf('\n\n');


fprintf('p value of a t-test between the sighted and LCA1 for each other measure for just those subjects with protocol 2 data:\n');
for i=1:3
    
        
tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexblind)=1;
indexblind_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexsight)=1;
indexsight_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexlca1)=1;
indexlca1_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);


    [~,p,~,stats]=ttest2(othermeasures(indexsight_protocol2,i),othermeasures(indexlca1_protocol2,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.3f') '\n']);
end
fprintf('\n\n');

fprintf('p value of a t-test between the blind and LCA1 for each other measure for just those subjects with protocol 2 data:\n');
for i=1:3
    
        
tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexblind)=1;
indexblind_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexsight)=1;
indexsight_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexlca1)=1;
indexlca1_protocol2=find(~isnan(othermeasures(:,i)).*(protocolNumber==2).*tmpSubjectArray);


    [~,p,~,stats]=ttest2(othermeasures(indexblind_protocol2,i),othermeasures(indexlca1_protocol2,i));
    fprintf(['measure ' LabelsMeasures{i} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.3f') '\n']);
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
tmpSubjectArray(indexblind)=1;
indexblind_protocol2=find(~isnan(faDataAgeRegressed).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexsight)=1;
indexsight_protocol2=find(~isnan(faDataAgeRegressed).*(protocolNumber==2).*tmpSubjectArray);

tmpSubjectArray=zeros(NumSubjects,1);
tmpSubjectArray(indexlca1)=1;
indexlca1_protocol2=find(~isnan(faDataAgeRegressed).*(protocolNumber==2).*tmpSubjectArray);



fprintf('p value of a t-test between the blind and sighted for each other measure:\n');
    [~,p,~,stats]=ttest2(faDataAgeRegressed(indexblind_protocol2),faDataAgeRegressed(indexsight_protocol2));
    fprintf(['measure ' LabelsMeasures{3} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.5f') '\n']);
fprintf('\n\n');


fprintf('p value of a t-test between the sighted and LCA1 for each other measure:\n');
    [~,p,~,stats]=ttest2(faDataAgeRegressed(indexsight_protocol2),faDataAgeRegressed(indexlca1_protocol2));
    fprintf(['measure ' LabelsMeasures{3} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.5f') '\n']);
fprintf('\n\n');

fprintf('p value of a t-test between the blind and LCA1 for each other measure:\n');
    [~,p,~,stats]=ttest2(faDataAgeRegressed(indexblind_protocol2),faDataAgeRegressed(indexlca1_protocol2));
    fprintf(['measure ' LabelsMeasures{3} ': ' 't-test (' num2str(stats.df) 'df): ' num2str(stats.tstat,'%1.2f') '  p=' num2str(p,'%1.5f') '\n']);
fprintf('\n\n');
