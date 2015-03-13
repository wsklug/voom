clear all;
close all;
clc;

% Read in nodal coordinates - all data on one line
NumFilesStart = 0;
NumFilesEnd = 1000;
DeltaNumFile = 10;

indCap = 0;
% indCapVST = 0;

for fnum = NumFilesStart:DeltaNumFile:NumFilesEnd
    clearvars -except fnum NumFilesStart DeltaNumFile NumFilesEnd indCap NumCap % CapVStemperature indCapVST
    NumCap(indCap+1,1:7) = 0;
    fnum
    Nodes = [];
    NBPent = 0; %8*2;
    % fid = fopen(['~/Desktop/Cyl13H/Cyl13Hiter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/U8sp/U8spiriter_Pr_',num2str(fnum),'.vtk'],'r');
    fid = fopen(['~/Desktop/Sphere/Sphiter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/Sphere/T7_6iter_Pr_',num2str(fnum),'.vtk'],'r');
%     fid = fopen(['~/Desktop/TestGeom/ParB/ParBiter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/FromUnduloids/MeltingT/U5spir/U5spirMeltingiter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/FromUnduloids/JLgeom/U8spirSub4iter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/FromUnduloids/U7sp/U7spirSub4iter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/FromUnduloids/MeltingT10B/JLunduloid4iter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/FromIchos/T13/T13_5iter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/FromIchos/T13/T13_5iter_Pr_',num2str(fnum),'.vtk'],'r');
    % fid = fopen(['~/Desktop/JLgeom/U8spirSub5iter_Pr_',num2str(fnum),'.vtk'],'r');
    header = fgetl(fid);
    header = fgetl(fid);
    header = fgetl(fid);
    header = fgetl(fid);
    header = fgetl(fid);
    NumNodes = str2num(header(8:12));
    
    % Read nodes
    for j = 1:NumNodes
        header = fgetl(fid);
        Nodes(j,:) = str2num(header);
    end
    
    % Add nodes at the center
    % Icshosahedral
%     Nodes(NumNodes+1,:) = mean(Nodes); % for T numbers
    % Unduloids and cylinders
%     Zmax = 1.0*(max(Nodes(:,3)));
%     Zmin = 1.0*(min(Nodes(:,3)));
%     Nadd = 100;
%     Nodes(NumNodes+1:NumNodes+Nadd,:) = [zeros(Nadd,1), zeros(Nadd,1), transpose(linspace(Zmin,Zmax,Nadd))];
    % Sphere
      Nodes(NumNodes+1,:) = [mean(Nodes)];
      %  Nodes(NumNodes+2:NumNodes*2+1,:) = [Nodes(1:NumNodes,:)*0.5]; No
      %  change in the posttprocessing results by adding more points in the
      %  sphere
    % Paraboloids
%     Ymax = 1.0*(max(Nodes(:,2)));
%     Ymin = 1.0*(min(Nodes(:,2)));
%     Nadd = 100;
%     Nodes(NumNodes+1:NumNodes+Nadd,:) = [zeros(Nadd,1), transpose(linspace(Ymin,Ymax,Nadd)), zeros(Nadd,1)];
    

    % Compute 3D Delaunnay triangulation 
    Tri = delaunay(Nodes(:,1),Nodes(:,2),Nodes(:,3));
%     tetramesh(Tri,Nodes);
    
    % Extract surface mesh
    NumEl = 0;
    SurfaceNodes = [1:1:NumNodes];
    for i=1:size(Tri,1)
        ElInt = intersect(Tri(i,:),SurfaceNodes);
        if (length(ElInt) == 3)
            NumEl = NumEl + 1;
            ConnTable(NumEl,:) = ElInt;
        end
    end
    
    
%     hold on;
%     for i = 1:size(ConnTable,1)
%             A = Nodes(ConnTable(i,1),:);
%             B = Nodes(ConnTable(i,2),:);
%             C = Nodes(ConnTable(i,3),:);
% 
%             plot3([A(1),B(1),C(1),A(1)],[A(2),B(2),C(2),A(2)],[A(3),B(3),C(3),A(3)],'k','Linewidth',1.5);
% 
%     end
%     hold off;
    

    
    % Write computed surface mesh - just to check the mesh is ok
    ElementsNumber = size(ConnTable,1);
%     % Write vtk file to imported into voom
%     % Printing vtk file
%     % Write vtk file
%     % fileID = fopen(['~/Desktop/FromUnduloids/MeltingT/U5spir/JL5spcapMesh_',num2str(fnum),'.vtk'],'w');
%     % fileID = fopen(['~/Desktop/FromUnduloids/JLgeom/JL8spcapMesh_',num2str(fnum),'.vtk'],'w');
%     % fileID = fopen(['~/Desktop/FromUnduloids/U7sp/U7spcapMesh_',num2str(fnum),'.vtk'],'w');
%     % fileID = fopen(['~/Desktop/FromUnduloids/MeltingT10B/JLcapMesh_',num2str(fnum),'.vtk'],'w');
%     % fileID = fopen(['~/Desktop/FromIchos/T13/T13capMesh_',num2str(fnum),'.vtk'],'w');
%     % fileID = fopen(['~/Desktop/FromIchos/T13/T13capMesh_',num2str(fnum),'.vtk'],'w');
%     % fileID = fopen(['~/Desktop/JLgeom/U8spirSubCapMesh_',num2str(fnum),'.vtk'],'w');
%     fprintf(fileID, '%s\n', '# vtk DataFile Version 3.0');
%     fprintf(fileID, '%s\n', 'vtk output');
%     fprintf(fileID, '%s\n', 'ASCII');
%     fprintf(fileID, '%s\n', 'DATASET POLYDATA');
%     fprintf(fileID, '%s ' , 'POINTS');
%     fprintf(fileID, '%d ' , NumNodes);
%     fprintf(fileID, '%s\n', 'float');
% 
%     for i = 1:NumNodes
%         fprintf(fileID, '%1.8f %1.8f %1.8f \n', Nodes(i,:));
%     end
% 
%     fprintf(fileID, '%s ' , 'POLYGONS');
%     fprintf(fileID, '%d ' , ElementsNumber);
%     fprintf(fileID, '%d\n', ElementsNumber*4);
% 
%     for i = 1:ElementsNumber
%         fprintf(fileID, '%d %d %d %d \n',[3, ConnTable(i,1:3)-1]);
%     end




    
    
    
    % For every vertex
    Outline = zeros(NumNodes,9);
    Cnum = 0;
    for i = 1:NumNodes
        % i
        Cnum = Cnum+1;
        ind = 0;
        Triangles = [];
        % Collect triangles sharing a vertex
        for j = 1:ElementsNumber
            if ( sum(ConnTable(j,:) == i) == 1)
                ind = ind + 1;
                Triangles(ind) = j;
            end
        end
        
%             hold on;
%             for ii = 1:length(Triangles)
%                     A = Nodes(ConnTable(Triangles(ii),1),:);
%                     B = Nodes(ConnTable(Triangles(ii),2),:);
%                     C = Nodes(ConnTable(Triangles(ii),3),:);
% 
%                     plot3([A(1),B(1),C(1),A(1)],[A(2),B(2),C(2),A(2)],[A(3),B(3),C(3),A(3)],'k','Linewidth',1.5);
% 
%             end
%             hold off;
        
        % Order triangles
        Tunder = ConnTable(Triangles(1),:);
        Attempt = 0;
        SwitchDir = 0;
        Outline(i,1) = Triangles(1);
        Triangles(1) = [];
        indOut = 2;
        while (length(Triangles) > 0)
            for j = 1:length(Triangles)
                a = intersect(Tunder, ConnTable(Triangles(j),:));
                Attempt = Attempt+1;
                if (length(a) == 2)
                    Tunder = ConnTable(Triangles(j),:);
                    Attempt = 0;
                    if (SwitchDir == 0)
                        Outline(i,indOut) = Triangles(j);
                    else
                        temp = Outline(i,1:indOut-1);
                        Outline(i,2:indOut) = temp;
                        Outline(i,1) = Triangles(j);
                    end
                    
                    Triangles(j) = [];
                    indOut = indOut +1;
                    Cnum = Cnum+1;
                    break;
                end
            end
%             b = intersect(ConnTable(Outline(i,1),:), ConnTable(Outline(i,indOut-1),:));
%             if (length(b) == 2 && length(Triangles) < 2)
%                 Triangles = [];
%             end
            if (length(Triangles) == 1)
                b = intersect(ConnTable(Outline(i,1),:), ConnTable(Triangles(1),:));
                if (length(b) == 2)
                    if (SwitchDir == 0)
                        Outline(i,indOut) = Triangles(1);
                    else
                        temp = Outline(i,1:indOut-1);
                        Outline(i,2:indOut) = temp;
                        Outline(i,1) = Triangles(1);
                    end
                    Triangles = [];
                    indOut = indOut +1;
                    Cnum = Cnum+1;
                    Attempt = 0;
                end
            end
            
            % If loop is not a complete loop and we did not start from an
            % exteme, next trianle may attach at the beginning rather than
            % at the end
            if (Attempt == length(Triangles))
                Tunder = ConnTable(Outline(i,1),:);
                SwitchDir = 1;
            end
            
        end
        
        Outline(i,indOut) = Outline(i,1);
        Cnum = Cnum+2;
        
    end
    
    % Center of triangles used for outlines
    for j = 1:ElementsNumber
        OutNodes(j,:) = mean(Nodes(ConnTable(j,:),:));
    end
    
 

    % Write outline file
    % fileID = fopen(['~/Desktop/Cyl13H/Cyl13Hcap_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/U8sp/U8spCap_',num2str(fnum),'.vtk'],'w');
    fileID = fopen(['~/Desktop/Sphere/SphereCap_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/Sphere/SphereCap_',num2str(fnum),'.vtk'],'w');
%     fileID = fopen(['~/Desktop/TestGeom/ParB/ParBcap_Pr_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/FromUnduloids/MeltingT/U5spir/MeltingU5spCap_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/FromUnduloids/JLgeom/JL8spCap_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/FromUnduloids/U7sp/U7spCap_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/FromUnduloids/MeltingT10B/JLgeomCap_',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/FromIchos/T13/T13cap',num2str(fnum),'.vtk'],'w');
    % fileID = fopen(['~/Desktop/FromIchos/T13/T13cap',num2str(fnum),'.vtk'],'w');   
    % fileID = fopen(['~/Desktop/JLgeom/U8spirSubCap_',num2str(fnum),'.vtk'],'w');
    fprintf(fileID, '%s\n', '# vtk DataFile Version 3.0');
    fprintf(fileID, '%s\n', 'vtk output');
    fprintf(fileID, '%s\n', 'ASCII');
    fprintf(fileID, '%s\n', 'DATASET POLYDATA');
    fprintf(fileID, '%s ' , 'POINTS');
    fprintf(fileID, '%d ' , size(OutNodes,1));
    fprintf(fileID, '%s\n', 'float');

    for i=1:size(OutNodes,1)
        fprintf(fileID, '%1.6f %1.6f %1.6f \n', OutNodes(i,:));
    end

    fprintf(fileID, '\n%s ' , 'POLYGONS');
    fprintf(fileID, '%d ' , size(Outline,1));
    fprintf(fileID, '%d\n', Cnum);

    indCap = indCap + 1;
    for i=1:size(Outline,1)
        if (Outline(i,4) == 0 ) 
            fprintf(fileID, '%d %d %d %d \n',[3, Outline(i,1:3)-1]);
            CellValue(i) = 2;
            NumCap(indCap, 1) = NumCap(indCap,1)+1;
        elseif (Outline(i,5) == 0 ) 
            fprintf(fileID, '%d %d %d %d %d \n',[4, Outline(i,1:4)-1]);
            CellValue(i) = 3;
            NumCap(indCap, 2) = NumCap(indCap,2)+1;
        elseif (Outline(i,6) == 0 ) 
            fprintf(fileID, '%d %d %d %d %d %d \n',[5, Outline(i,1:5)-1]);
            CellValue(i) = 4;
            NumCap(indCap, 3) = NumCap(indCap,3)+1;
        elseif (Outline(i,7) == 0 ) % pentamers
            fprintf(fileID, '%d %d %d %d %d %d %d \n',[6, Outline(i,1:6)-1]);
            CellValue(i) = 5;
            NumCap(indCap, 4) = NumCap(indCap,4)+1;
        elseif (Outline(i,8) == 0 ) % hexamers
            fprintf(fileID, '%d %d %d %d %d %d %d %d \n',[7, Outline(i,1:7)-1]);
            CellValue(i) = 6;
            NumCap(indCap, 5) = NumCap(indCap,5)+1;
        elseif (Outline(i,9) == 0 ) % heptamers
            fprintf(fileID, '%d %d %d %d %d %d %d %d %d \n',[8, Outline(i,1:8)-1]);
            CellValue(i) = 7;
            NumCap(indCap, 6) = NumCap(indCap,6)+1;
        else  
            fprintf(fileID, '%d %d %d %d %d %d %d %d %d %d \n',[9, Outline(i,1:9)-1,]);
            CellValue(i) = 8;
            NumCap(indCap, 7) = NumCap(indCap,7)+1;
        end
    end
   
    fprintf(fileID, '\n%s ' , 'CELL_DATA');
    fprintf(fileID, '%d \n' , size(Outline,1));
    fprintf(fileID, '%s \n' , 'SCALARS    Valence    double    1');
    fprintf(fileID, '%s \n' , 'LOOKUP_TABLE default');

    for i=1:size(Outline,1)
        fprintf(fileID, '%d \n',CellValue(i));
    end
    
    NumCap(indCap, :)
    
%     if(indCap == 50)
%         indCapVST = indCapVST+1;
%         CapVStemperature(indCapVST, :) = mean(Numcap);
%         indCap = 0;
%     end

    fclose(fileID);
    
end



% Compute average mesh size in last step
rm = 0.0;
NumEl = size(ConnTable,1);
for i=1:NumEl
    Element = ConnTable(i,:);
    rm = rm + ( norm(Nodes(Element(2),:)-Nodes(Element(1),:)) + ....
                norm(Nodes(Element(3),:)-Nodes(Element(2),:)) + ....
                norm(Nodes(Element(1),:)-Nodes(Element(3),:)) )/3.0;
end
rm = rm/NumEl
sigma = rm/(2^(1/6))

figure(1)
hold on;
set(gca,'FontSize',14)
StepConstT = 100;
plot([StepConstT],NumCap(2,4), '-sb','Linewidth', 2, 'MarkerSize', 5);
plot([StepConstT],NumCap(2,5), '-ok','Linewidth', 2, 'MarkerSize', 5);
plot([StepConstT],NumCap(2,6), '-vr','Linewidth', 2, 'MarkerSize', 5);

NumCapSize = size(NumCap,1);
plot([0:StepConstT:StepConstT*(NumCapSize-1)], NumCap(:,4), '-b');
plot([0:StepConstT:StepConstT*(NumCapSize-1)], NumCap(:,5), '-k');
plot([0:StepConstT:StepConstT*(NumCapSize-1)], NumCap(:,6), '-r');

plot([0:StepConstT*10:(NumCapSize-1)*StepConstT], NumCap(1:10:NumCapSize,4), 'sb', 'Linewidth', 3, 'MarkerSize', 5);
plot([0:StepConstT*10:(NumCapSize-1)*StepConstT], NumCap(1:10:NumCapSize,5), 'ok', 'Linewidth', 3, 'MarkerSize', 5);
plot([0:StepConstT*10:(NumCapSize-1)*StepConstT], NumCap(1:10:NumCapSize,6), 'vr', 'Linewidth', 3, 'MarkerSize', 5);
xlabel('MC steps');
ylabel('Capsomers');
legend('Pentamers','Hexamers','Heptamers');
xlim([0 10000]);
hold off;
    
figure(2)
hold on;
set(gca,'FontSize',14)
NumCapSize = size(NumCap,1);
plot([1:size(NumCap,1)], NumCap(:,4)-NumCap(:,6)-NBPent, '-k');
plot([1:size(NumCap,1)], NumCap(:,4)-NBPent, '-b');
plot([1:size(NumCap,1)], NumCap(:,6), '-r');
StepConstT = 10;
plot([StepConstT:StepConstT:NumCapSize], NumCap(StepConstT:StepConstT:NumCapSize,4)-NumCap(StepConstT:StepConstT:NumCapSize,6)-NBPent, 'sk', 'Linewidth', 3, 'MarkerSize', 5);
plot([StepConstT:StepConstT:NumCapSize], NumCap(StepConstT:StepConstT:NumCapSize,4)-NBPent, 'sb', 'Linewidth', 3, 'MarkerSize', 5);
plot([StepConstT:StepConstT:NumCapSize], NumCap(StepConstT:StepConstT:NumCapSize,6), 'sr', 'Linewidth', 3, 'MarkerSize', 5);
xlabel('MC steps');
ylabel('Capsomers');
legend('Pentamers - Heptamers','Pentamers','Heptamers');
xlim([0 100]);
hold off;
        
% save('~/Desktop/Cyl13H/CapsAnnCyl13H.dat','NumCap','-ascii');
% save('~/Desktop/U8sp/CapsAnnU8sp.dat','NumCap','-ascii');
save('~/Desktop/Sphere/CapsAnnSphere.dat','NumCap','-ascii');




