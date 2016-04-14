clear all;
close all;
clc;


% Parameters to run the code
Periodic      = 0;
NumFilesStart = 0;
NumFilesEnd   = 1500;
DeltaNumFile  = 1;
Length        = 7.74094;
dirname       = '/Users/luigi/Documents/Extra/Archaea/Harmonic/H';


%% Construct conn table for spring particles
SpiralCenters = [912, 916, 914, 918, 913, 915, 917, 919];

ind = 1;
for j = 0:7 % 8 spirals
  sh = j*114;
  for i = 0:55 
    HarmonicConn(ind, 1) = i + sh;
    HarmonicConn(ind, 2) = i+1 + sh;
    ind = ind + 1;
  end
  for i = 0:55
    HarmonicConn(ind, 1) = 113 - i + sh;
    HarmonicConn(ind, 2) = 113 - i - 1 + sh;
    ind = ind + 1;
  end
  
  HarmonicConn(ind, 1) = 56 + sh;
  HarmonicConn(ind, 2) = SpiralCenters(j+1);
  ind = ind + 1;

  HarmonicConn(ind, 1) = SpiralCenters(j+1);
  HarmonicConn(ind, 2) = 113 + sh;
  ind = ind+1;
end



%% Read in nodal coordinates - all data on one line


for folder = 2:3
    indCap  = 0;
    indFile = 1;    
    for fnum = NumFilesStart:DeltaNumFile:NumFilesEnd
        fnum
        Nodes = [];

        fid = fopen([dirname,num2str(folder),'/U8spiriter_Pr_',num2str(fnum),'.vtk'],'r');
        
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

        fclose(fid);
    
       
 
        
        %% Print Endo surface to vtk file
        fileID = fopen([dirname,num2str(folder),'/HarmSpring_',num2str(fnum),'.vtk'],'w');

        fprintf(fileID, '%s\n', '# vtk DataFile Version 3.1');
        fprintf(fileID, '%s\n', 'vtk output');
        fprintf(fileID, '%s\n', 'ASCII');
        fprintf(fileID, '%s\n', 'DATASET UNSTRUCTURED_GRID');
        fprintf(fileID, '%s ' , 'POINTS');
        fprintf(fileID, '%d ' , size(Nodes,1));
        fprintf(fileID, '%s\n', 'float');

        for i = 1:size(Nodes,1)
            fprintf(fileID, '%1.6f %1.6f %1.6f \n', Nodes(i,1:3));
        end

        fprintf(fileID, '\n%s ' , 'CELLS');
        fprintf(fileID, '%d ' , size(HarmonicConn,1));
        fprintf(fileID, '%d\n', (1 + size(HarmonicConn,2))*size(HarmonicConn,1) );

        for i=1:size(HarmonicConn,1)
            fprintf(fileID, '%d %d %d %d \n',[2, HarmonicConn(i,1:2)]);
        end   

        fprintf(fileID, '\n%s ' , 'CELL_TYPES');
        fprintf(fileID, '%d \n' , size(HarmonicConn,1));
        for i=1:size(HarmonicConn,1)
            fprintf(fileID, '%d\n', 3);
        end    

        fclose(fileID);
        
    end
    
end



