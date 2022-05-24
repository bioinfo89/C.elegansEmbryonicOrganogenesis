function script31_predict_cell_position_for_figure(env,prjRoot,dataMode,testMode,axisOri,plotStartNumCell,plotEndNumCell,trainTStart,trainTEnd,outputPathVideo,NumTrees)
close all;

countFrame = 1;
F(1).cdata=[];
F(1).colormap=[];

corMarkersize=6

for t = trainTStart:trainTEnd

    trainFN = dir(['D:\home\bxn200003\scratch\test_zone\avg_embryo\trainData_shuffled\train_time_' num2str(t) '_numCell_*.csv']);
    if size(trainFN,1)>0
        trainFN =['D:\home\bxn200003\scratch\test_zone\avg_embryo\trainData_shuffled\' trainFN.name];
        
        ds = importdata(trainFN)

        cellNames = ds.textdata(:,1);
        cellNames = cellNames(2:end,1);

        coord=ds.data(:,end-2:end);
        geneExpr = ds.data(:,1:end-3);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Do random forest regression along Y given gene express as X
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        numCell = size(cellNames,1);
         
        X=[]
        Y=[];
        if strcmp(dataMode,'gene')
            X=geneExpr(:,:);         
            if strcmp(axisOri,'A-P')
                Y=coord(:,2);
            end
            if strcmp(axisOri,'L-R')
                Y=coord(:,1);
            end
            if strcmp(axisOri,'V-D')
                Y=coord(:,3);
            end
        else
            X=rand(size(geneExpr,1),size(geneExpr,2));
            if strcmp(axisOri,'A-P')
                Y=coord(:,2);
            end
            if strcmp(axisOri,'L-R')
                Y=coord(:,1);
            end
            if strcmp(axisOri,'V-D')
                Y=coord(:,3);
            end
        end      
        
        Xtrain = X(1:2:end,:);
        Ytrain=Y(1:2:end,:);
        
        Xtest = X(2:2:end,:);
        Ytest=Y(2:2:end,:);
        
        XcleanTrain=[];        
        for p=1:size(Xtrain,2)
            v=Xtrain(:,p);
            if sum(v==-100)~=size(v,1)
                XcleanTrain=[XcleanTrain,v];
            end            
        end
          
        XcleanTest=[];
        for p=1:size(Xtest,2)
            v=Xtest(:,p);
            if sum(v==-100)~=size(v,1)
                XcleanTest=[XcleanTest,v];
            end            
        end
        
        h1=figure(1)
        set(h1,'position',[680 74 1213 904])       
        if numCell<=plotEndNumCell & numCell>plotStartNumCell
            % do regression
            
            B = TreeBagger(NumTrees,XcleanTrain,Ytrain,'method','regression')
            
            predX=[];
            if strcmp(testMode,'train')
               predX=XcleanTrain;
            else                
               predX=XcleanTest;
            end
            Y
            yPredict = predict(B,predX)

            subplot(2,2,[1,2])
            if strcmp(testMode,'test')
                %plot(Ytest,yPredict,'ko'); hold off
                
                minYtest = min(Ytest);
                maxYtest = max(Ytest);
                minyPredict = min(yPredict);
                maxyPredict = max(yPredict);
                
                set(gcf,'position',[50   555   560   420])
                Z = getKDE2D(Ytest',yPredict',minYtest,maxYtest,minyPredict,maxyPredict,8);

                imagesc(Z)
                hold on
                colormap summer
                contour(Z,3,'LineColor','w','color',[0.7 0.7 0.7],'linewidth',2);
                axis square

                hold on

                plot(300*(Ytest-minYtest)/(maxYtest-minYtest),200*(yPredict-minyPredict)/(maxyPredict-minyPredict),'wo','markerfacecolor','r','linewidth',1,'markersize',corMarkersize,'color','w'); hold on

                dim1 = 300*(Ytest-minYtest)/(maxYtest-minYtest);
                dim2 = 200*(yPredict-minyPredict)/(maxyPredict-minyPredict);

                corel = corrcoef(dim1,dim2)
                text(15,188,['Pearson corr.= ' num2str(ceil((1000*corel(1,2)))/1000)],'color','w', 'fontsize',11,'fontweight','bold','FontName', 'Arial')
            
                hold off
                set(gca,'Ydir','normal')
                set(gca,'xtick',[])
                set(gca,'ytick',[])
                xlabel(['Observed cell position (' axisOri ' axis)'],'fontsize',12,'fontweight','bold');
                ylabel(['Predicted cell position'],'fontsize',12,'fontweight','bold');
                title(['T = ' num2str(t)  ', Num. of cell = ' num2str(numCell)],'fontsize',14);
                set(gcf,'color','w')
                axis square        
               
            else
                plot(Ytrain,yPredict,'ko'); hold off
            end
            
            xlabel(['Observed cell position along ' axisOri ' axis'],'fontsize',12);
            ylabel(['Predicted cell position along ' axisOri ' axis'],'fontsize',12);
            title(['T = ' num2str(t)  ', Num. of cell = ' num2str(numCell)],'fontsize',14);
            set(gcf,'color','w')
            colorbar 
            axis square
            pause(0.1)
            
            h1000 = figure(1000)
            set(gcf,'position',[50   658   403   317])
            
            minYtest = min(Ytest);
            maxYtest = max(Ytest);
            minyPredict = min(yPredict);
            maxyPredict = max(yPredict);
                
            set(gcf,'position',[50   555   560   420])
            Z = getKDE2D(Ytest',yPredict',minYtest,maxYtest,minyPredict,maxyPredict,8);

            imagesc(Z)
            hold on
            colormap summer
            contour(Z,3,'LineColor','w','color',[0.7 0.7 0.7],'linewidth',2);
            axis square

            hold on

            plot(300*(Ytest-minYtest)/(maxYtest-minYtest),200*(yPredict-minyPredict)/(maxyPredict-minyPredict),'wo','markerfacecolor','r','linewidth',1,'markersize',corMarkersize,'color','w'); hold on

            dim1 = 300*(Ytest-minYtest)/(maxYtest-minYtest);
            dim2 = 200*(yPredict-minyPredict)/(maxyPredict-minyPredict);
            
            corel = corrcoef(dim1,dim2)
            text(15,188,['Pearson corr.= ' num2str(ceil((1000*corel(1,2)))/1000)],'color','w', 'fontsize',11,'fontweight','bold','FontName', 'Arial')
            
            hold off
            set(gca,'Ydir','normal')
            ax=gca;
            set(gca,'xtick',[])
            set(gca,'ytick',[])
            ax.XLabel.FontSize = 12;
            ax.YLabel.FontSize = 12;
            xlabel(['Observed cell position (' axisOri ' axis)'],'fontsize',10,'fontweight','bold');
            if strcmp(dataMode,'rand')
                ylabel(['Predicted cell position (Neg. ctrl.)'],'fontsize',11,'fontweight','bold');
            else
                ylabel(['Predicted cell position'],'fontsize',11,'fontweight','bold');                
            end
            title(['T = ' num2str(t)  ', Num. of cell = ' num2str(numCell)],'fontsize',12);
            set(gcf,'color','w')
            colorbar
            axis square        
            set(gcf,'position',[50   645   363   330])
     
            h1=figure(1)
            
            saveas(h1000,[outputPathVideo 'reference_embryo_dataMode_' dataMode '_testMode_' testMode '_axisOri_' axisOri '_plotStartNumCell_' num2str(plotStartNumCell) '_plotEndNumCell_' num2str(plotEndNumCell) '_time_' num2str(t) '.fig']);
            saveas(h1,[outputPathVideo 'reference_embryo_dataMode_' dataMode '_testMode_' testMode '_axisOri_' axisOri '_plotStartNumCell_' num2str(plotStartNumCell) '_plotEndNumCell_' num2str(plotEndNumCell) '_time_' num2str(t) '_h1.fig']);
           
            save([outputPathVideo 'reference_embryo_dataMode_' dataMode '_testMode_' testMode '_axisOri_' axisOri '_plotStartNumCell_' num2str(plotStartNumCell) '_plotEndNumCell_' num2str(plotEndNumCell) '_time_' num2str(t) '.mat'], 'corel','dim1','dim2');

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
	  set(h1,'position',[680 74 1213 904])
        ax1=subplot(2,2,3)
        plot3(coord(:,1),coord(:,2),coord(:,3),'ro');
        for i=1:size(coord,1)
            text(coord(i,1),coord(i,2),coord(i,3),char(cellNames{i})); hold on
        end

        xlabel('V-D axis');ylabel('A-P axis');zlabel('L-R axis');
        
        title(['Reference Embryo: T=' num2str(t) ', Number cell=' num2str(numCell)]);
        grid on
        ax = gca;
       
        set(gca,'DataAspectRatio',[1 1 1])
        axis([-17 17 -30 30 -17 17])
        hold off
        
        ax2=subplot(2,2,4)
        cla(ax2) 

        numcell=size(coord,1)
        c = coord;
        
        a = (coord(:,2)-min(coord(:,2)))./(max(coord(:,2))-min(coord(:,2)));
        a = a*0.5+0.5;
        color = [zeros(size(a,1),1),a,zeros(size(a,1),1)];
        
        r=ones(1,numcell)*1.5
        alpha = 1

        bubbleplot3(c(:,1),c(:,2),c(:,3),r,color,alpha);
        
        camlight left; view(60,30);
        camlight('headlight') ; view(60,30);
        camlight('left') ; view(60,30);
        lighting gouraud
        
        light('Position',[-30 -10 -10],'Style','local')
        light('Position',[-30 10 -10],'Style','local')
        light('Position',[-30 0 0],'Style','local')
        
        grid off
   
        ax = gca;
        %ax.BoxStyle = 'full';
        box off
        pause(1)
        
        set(gca,'DataAspectRatio',[1 1 1])
        xlabel('V-D axis','fontweight','bold','fontsize',14);ylabel('A-P axis','fontweight','bold','fontsize',14);zlabel('L-R axis','fontweight','bold','fontsize',14);
        
        title(['Reference Embryo: T=' num2str(t)  ', Number cell=' num2str(numCell)],'fontweight','bold','fontsize',16);        
       
        linkprop([ax2 ax1], {'View', 'XLim', 'YLim', 'ZLim',})
        
        hold off
        
        set(gcf,'color','w')
        axis([-17 17 -30 30 -17 17])
        
        hold off
        pause(0.5)
        
        F(countFrame) = getframe(gcf);             
        countFrame = countFrame+1;     
        end 
    end
end

% write the frames to the video
writerObj = VideoWriter([outputPathVideo 'reference_embryo_dataMode_' dataMode '_testMode_' testMode '_axisOri_' axisOri '_plotStartNumCell_' num2str(plotStartNumCell) '_plotEndNumCell_' num2str(plotEndNumCell) '.avi']);
writerObj.FrameRate = 2;
open(writerObj);               
for v=1:length(F)
   frame = F(v) ;    
   if size(frame.cdata,1)>0
    writeVideo(writerObj, frame);
   end
end 
close(writerObj);  
clear F


