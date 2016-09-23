%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script will give rotation matrices to rotate the T7input.vtk viral
%capsid to orient its three-fold and two-fold axes along the Z axes so that
%our indentation code can indent along each of these axes without modifying
%the code itself.
%
%Author: Amit Singh
%Year: 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%Read all the (x,y,z) co-ordinates into 'allPoints'
allPoints = dlmread('allPoints.dat','\t');

%Input Valence matrix obtained from insertVtkInValence() in voom
fin = fopen('valence.dat','r');
fout = fopen('valence2.dat','w');
while ~feof(fin)
        s = fgetl(fin);
        s = strrep(s,'\n','\t');
        fprintf(fout,'%s',s);
end

fclose(fin);
fclose(fout);
movefile('valence2.dat','valence.dat');

valence =dlmread('valence.dat','\t');

%Now we will use 'valence' vector to identify the node numbers that
%correspond to pentamers and then we will extract their (x,y,z)
%co-ordinates from 'allPoints' matrix and store it in 'penta'.
[~,pentaRows] = find(valence == 5);

%Read pentamer co-ordinates into 'penta' from 'allPoints'.
penta = allPoints(pentaRows,:);

%We arbitrarily choose one pentamer as reference for further calculations.
%In this case we will choose the one that has highest z-co-ordinate.
%Variable refRow contains the row number of the pentamer with highest z.
[refRow,~] = find(penta(:,3) == max(penta(:,3)));
refRow = refRow(1);

%Now we want to calculate the Eucledian distance of all other pentamers
%with respect to the reference pentamer identified above. We expect to find
%5 equidistant pentmers that are closest to the reference.
%bsxfun() is a good way to avoid writing for-loops over matrix elements
%that I found today

%First calculate a 12x3 matrix whose columns are (x-xref)^2, (y-yref)^2 and
%(z-zref)^2 where (xref,yref,zref) are co-ordinates of the reference
%pentamer
diffSquares = bsxfun(@power,bsxfun(@minus,penta,penta(refRow,:)),2);

%The next line returns a 12x1 matrix of sum of squares
%(x-xref)^2+(y-yref)^2+(z-zref)^2
sumDiffSquares = sum(diffSquares,2);

%Finally, the elementwise square-root to give distances
distances = sqrt(sumDiffSquares);

%We can use the norm() to calculate sqrt((x-xref)^2+(y-yref)^2+(z-zref)^2)


%Now we have to find the row numbers from 'distances' of the pentamers that
%are closest to the reference. Of course, for row number equal to 'refRow'
%the distance is 0 because it is the reference pentamer itself. So we have
%to find the minimum value that is greater than 0. We expect 5 of them.

%We will replace the value distance(refRow,1) with some large value so that
%we can use min() to identify the minimum distances
distances(refRow,1) = sum(distances);
tol = 1e-5; %Tolerance
[closestRows,~] = find(abs(distances - min(distances)) <= tol);

%Consistency check: we should have only 5 closest pentamers
if (length(closestRows) ~= 5)
    error('Reference pentamer has %d nearest pentamers instead of 5.\n',...
        length(closestRows));
end

%Now we need to make 5 vectors going from the reference pentamer to each of
%the 5 closest pentamer.
vectors = bsxfun(@minus,penta(closestRows,:),penta(refRow,:));

%Now we need to choose any one (we will choose the first one) of the 5
%vectors we created above and find at least one of the remaining 4 vectors
%that is at an angle approximately 60degrees from it.
%There are 2 such vectors but finding 1 is good enough for us.
chosenVec = vectors(1,:);
dotProducts = vectors*chosenVec';
norms = sqrt(sum(vectors.^2,2));
normProds = norms*norm(chosenVec);
cosTheta = dotProducts./normProds;
thetaInDegrees = acos(cosTheta)*180/pi;

[row60deg,~] = find(abs(thetaInDegrees - 60) < tol);
thirdPentamerRow = row60deg(1,1);

%Three pentamers that form an equilateral face of icosahedron are given by
%the reference pentamer, the pentamer corresponding to the row index in
%'penta' whose value is given by first element of 'closestRows'. The value
%of 'thirdPentamer' gives the row number of 'closestRows' where we can find
%the row index for third pentamer from 'penta'
fP = penta(refRow,:); %First Pentamer coordinates
sP = penta(closestRows(1,1),:); %Seccond Pentamer coordinates
tP = penta(closestRows(thirdPentamerRow,1),:); %Third Pentame coordinates

%%%%%%%%%%%%%%%%%%% Rodrigues Rotation Matrix Formula %%%%%%%%%%%%%%%%%%%%%

% v_rot = R*v where v is the original vector and v_rot is the rotated
% vector. R is the required rotation matrix.
%
% R = I + sin(theta)K + (1-cos(theta))K^2 where
%
% K = [0,-k3,k2; k3,0,-k1; -k2,k1,0] where (k1,k2,k3) is the unit vector of
% the axis of rotation and theta is the angle of rotation
%
% For us, v=[0,0,1] the z-axis
% We need to calculate v_rot for the three-fold and two-fold axes.

%########################## For Three-fold axis ###########################
% Calculate v_rot for three-fold axes:
% We know it is perpendicular to the triangle formed by the three pentamers
% identified above and passes through the centroid of the triangle. It
% should point away from the origin.

triCentroid = (fP+sP+tP)/3;  %Surface normal to triangle through centroid
triCentroid = triCentroid/norm(triCentroid);

% Consistency check: Is the surface normal orthogonal to sides of the
% triangle
if (abs(acos(dot(triCentroid,fP-sP)/...
        (norm(triCentroid)*norm(fP-sP)))*180/pi - 90) > tol)
    error('Calculated three fold axis is not normal to triangular' +...
        'face of icosahedron!\n');
end
if (abs(acos(dot(triCentroid,fP-tP)/...
        (norm(triCentroid)*norm(fP-tP)))*180/pi - 90) > tol)
    error('Calculated three fold axis is not normal to triangular' +...
        'face of icosahedron!\n');
end

%Consistency check: Let's hardcode triCentroid as [0,1,0] i.e. y-axis. Do
%we get k=[1,0,0] and theta = 90 degrees? Just uncomment the next line and
%print k, theta and R3 at end of execution of script

%triCentroid = [0,1,0];

% Find k
v = [0,0,1];

k = cross(triCentroid,v);
k = k/norm(k);
K = [0,-k(3),k(2); k(3),0,-k(1); -k(2),k(1),0];

% Find theta
cosTheta = dot(v,triCentroid); % Both v and triCentroid are unit vectors
theta = acos(cosTheta);

I = [1,0,0;0,1,0;0,0,1]; %Identity matrix

%Finally, the rotation matrix for transforming v_rot to v. Note: we have
%used -theta instead of theta because theta takes v to v_rot and we want
%the other way round.
R3 = I + sin(-theta)*K + (1-cos(-theta))*K^2;

%Consistency check: R3 should belong to SO(3) i.e. det(R3)=1 and
%R3*transpose(R3)=I
if(abs(det(R3)-1)>tol || norm(R3*R3'-I)>tol)
    error('R3 does not belong to SO(3)!\n');
end

%########################## For Two-fold axis ############################
% Calculate v_rot for two-fold axis. We just want the midpoint of any side
% of the triangular face.
midPoint = fP + (sP-fP)/2;

% Consistency check: Is the midPoint position vector normal to side of the
% triangle
if (abs(acos(dot(midPoint,sP-fP)/...
        (norm(midPoint)*norm(sP-fP)))*180/pi - 90) > tol)
    error('Calculated two fold axis is not normal to side' +...
        'of the face of icosahedron!\n');
end

mP = midPoint/norm(midPoint);

% Find k
k = cross(mP,v);
k = k/norm(k);
K = [0,-k(3),k(2); k(3),0,-k(1); -k(2),k(1),0];

% Find theta
cosTheta = dot(v,mP); % Both v and triCentroid are unit vectors
theta = acos(cosTheta);

%Finally, the rotation matrix for transforming v_rot to v. Note: we have
%used -theta instead of theta because theta takes v to v_rot and we want
%the other way round.
R2 = I + sin(-theta)*K + (1-cos(-theta))*K^2;

%Consistency check: R2 should belong to SO(3) i.e. det(R2)=1 and
%R2*transpose(R2)=Id
if(abs(det(R2)-1)>tol || norm(R2*R2'-I)>tol)
    error('R2 does not belong to SO(3)!\n');
end

%########################## For Five-fold axis ############################
% The centre of any pentamer gives the position vector for five-fold axis
% We will choose the first one that we selected using refRows
fP = fP/norm(fP);
k = cross(fP,v);
k = k/norm(k);
K = [0,-k(3),k(2); k(3),0,-k(1); -k(2),k(1),0];

% Find theta
cosTheta = dot(v,fP); % Both v and triCentroid are unit vectors
theta = acos(cosTheta);

%Finally, the rotation matrix for transforming v_rot to v. Note: we have
%used -theta instead of theta because theta takes v to v_rot and we want
%the other way round.
R5 = I + sin(-theta)*K + (1-cos(-theta))*K^2;

%Consistency check: R5 should belong to SO(3) i.e. det(R5)=1 and
%R5*transpose(R5)=I
if(abs(det(R5)-1)>tol || norm(R5*R5'-I)>tol)
    error('R5 does not belong to SO(3)!\n');
end


%Now we will multiply them by R2, R3 and R5 to get the rotated points such
%that two-fold axis and three-fold axis are along z-direction respectively
twoFoldZ = allPoints*R2;
twoFold = fopen('twoFold.dat','w');
for l=1:size(twoFoldZ,1)
    fprintf(twoFold,'%18.16f\t%18.16f\t%18.16f\n',twoFoldZ(l,1),...
        twoFoldZ(l,2),twoFoldZ(l,3));
end
fclose(twoFold);

threeFoldZ = allPoints*R3;
threeFold = fopen('threeFold.dat','w');
for l=1:size(threeFoldZ,1)
    fprintf(threeFold,'%18.16f\t%18.16f\t%18.16f\n',threeFoldZ(l,1),...
        threeFoldZ(l,2),threeFoldZ(l,3));
end
fclose(threeFold);

fiveFoldZ = allPoints*R5;
fiveFold = fopen('fiveFold.dat','w');
for l=1:size(fiveFoldZ,1)
    fprintf(fiveFold,'%18.16f\t%18.16f\t%18.16f\n',fiveFoldZ(l,1),...
        fiveFoldZ(l,2),fiveFoldZ(l,3));
end
fclose(fiveFold);

