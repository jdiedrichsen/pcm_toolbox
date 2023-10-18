function varargout=pcm_recipe_correlation_MDS
% This example deomstrates the use of an MDS plot for the movements of 5 fingers for the left and
% right hand, respectively. 
% 
% 


load data_recipe_correlation.mat

% --------------------------------------
% 1. Crossvalidated covariance martric and correlation
for p=1:12
    Z=pcm_indicatorMatrix('identity',condVec{p});
    % Subtract mean for each hand and run
    Gcv(:,:,p)=pcm_estGCrossval(Data{p},partVec{p},condVec{p});
end;
mGcv = mean(Gcv,3); 

C1=pcm_indicatorMatrix('allpairs',[1:10]); 
c=pcm_indicatorMatrix('allpairs',[1:5]); 
C2=[c c];
C3= [eye(5) zeros(5,5)]; 
[Y,l]=pcm_classicalMDS(mGcv,'contrast',C1');

color={'r','b'};
for h=[1:2]
    indx = [1:5]+(h-1)*5; 
    plot3(Y(indx,1),Y(indx,2),Y(indx,3),'o','Color',color{h},'MarkerSize',10,'LineWidth',2); 
    hold on; 
    D = Y(indx,1:3); 
    D(end+1,:)=D(1,:); 
    line(D(:,1),D(:,2),D(:,3),'Color',color{h},'LineWidth',2); 
end; 
plot3(0,0,0,'k+','MarkerSize',10,'LineWidth',2); 
set(gca,'XGrid',1,'YGrid',1,'ZGrid',1,'XTickLabel',[],'YTickLabel',[],'ZTickLabel',[]);
hold off; 
set(gcf,'PaperPosition',[2 2 6 5])
axis equal;
wysiwyg; 
keyboard; 