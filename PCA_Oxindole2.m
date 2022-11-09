% Program to Perform PCA
clear
%load CMatrix.txt; % C matrix on calibration set without data rotation
load CMatrix_62.txt; % C matrix on calibration set without data rotation
load RMatrix1814_2.txt; % R matrix on calibration set without data rotation 1st cross exam
load CMatrix1814_2.txt; % C matrix on calibration set without data rotation 1st cross exam
%load CMatrix1814_2rot30.txt; % C matrix on calibration set with 30 degrees data rotation 1st cross exam
%load CMatrix1814_2rot60.txt; % C matrix on calibration set with 60 degrees data rotation 1st cross exam
%load CMatrix1814_2rot90.txt; % C matrix on calibration set with 90 degrees data rotation 1st cross exam
%load CMatrix1814_2ndCross.txt; % C matrix on calibration set without data rotation 2nd cross exam
%load CMatrix1814_3rdCross.txt; % C matrix on calibration set without data rotation 2nd cross exam
%load CMatrix1814_4thCross.txt; % C matrix on calibration set without data rotation 2nd cross exam
%load CMatrixrot.txt; % 50 degrees orthogonal
%load OxTrSet.txt; load OxTrSet2.txt; load OxTrSet3.txt; load OxTrSet4.txt; load OxTrSet5.txt;
%load OxTrSet6.txt;
load OxTrSet7.txt;
D1=OxTrSet7;
[m1,n1]=size(D1);

for i=1:m1
    for j=1:n1
        %D2(i,j)=(D1(i,j)-mean(D1(i,:)))/std(D1(i,:));
        D2(i,j)=(D1(i,j)-mean(D1(:,j)))/std(D1(:,j));
    end
end

% for OxTrSet4
%Dhl=D2(:,1:6);
%Dpol=D2(:,7:12);
%Dchar=D2(:,13:22);
%Dfukuip=D2(:,23:32);
%Dfukuin=D2(:,33:42);
%Dhf=D2(:,43);
%Dzpe=D2(:,44);
%Dtheren=D2(:,45);
%Dcv=D2(:,46);
%Dent=D2(:,47);
%Drot=D2(:,48:50);
%Dmm=D2(:,51);

% for OxTrSet6
Dhl=D2(:,1:6);
Dpol=D2(:,7:12);
Dhf=D2(:,13);
Dzpe=D2(:,14);
Dtheren=D2(:,15);
Dcv=D2(:,16);
Dent=D2(:,17);
Drot=D2(:,18:20);
Dmm=D2(:,21);

% for OxTrSet4
%D3(:,1:6)=Dhl; %D3(:,7:12)=Dpol; D3(:,13)=Dhf; D3(:,14)=Dzpe; D3(:,15)=Dtheren; 
%D3(:,1:6)=Dhl; D3(:,7)=Dhf; D3(:,8)=Dzpe; D3(:,9)=Dtheren; D3(:,10)=Dcv; D3(:,11)=Dent; D3(:,12:14)=Drot; %D3(:,15)=Dmm;
%D3(:,1:6)=D2(:,1:6); D3(:,7:14)=D2(:,43:50); D3(:,15)=Dmm; %D3(:,16:21)=Dpol;

% for OxTrSet6
%D3=D2;
D3(:,1:6)=Dhl; D3(:,7)=Dhf; D3(:,8)=Dzpe; D3(:,9)=Dtheren; D3(:,10)=Dcv; D3(:,11)=Dent; D3(:,12:14)=Drot; D3(:,15)=Dmm; %D3(:,16:21)=Dpol;

% remove outliers
%D4(1:21,:)=D3(2:22,:); D4(22:27,:)=D3(24:29,:); D4(28:50,:)=D3(31:53,:); 
%D4(51:52,:)=D3(55:56,:); D4(53:55,:)=D3(58:60,:);

% Training Set (row number in D3)
% cross examination study 1
D4=D3([2,4,6,9,11,14,16,17,18,22,24,27,28,31,34,37,39,42,43,45,46,48,49,50,53,54,56,57,58,60,61,62],:);
% anticancer (3, 7, 8, 10, 12, 13, 15, 20, 21, 25, 26, 29, 32, 33, 35, 36, 38, 40, 41) oxtrset7
% antimicrobial (44, 47, 51, 52, 55, 59, 63) oxtrset7

% cross examination study 2
%D4=D3([3,4,8,10,11,13,15,17,20,21,24,26,27,32,33,36,40,42,44,45,46,47,50,52,53,55,56,57,58,59,61,63],:);
% anticancer (2, 6, 7, 9, 12, 14, 16, 18, 22, 25, 28, 29, 31, 34, 35, 37, 38, 39, 41) oxtrset7
% antimicrobial (43, 48, 49, 51, 54, 60, 62) oxtrset7

% cross examination study 3
%D4=D3([2,6,7,9,12,13,14,16,20,21,25,26,28,32,34,36,38,41,43,44,46,47,49,50,51,53,54,56,57,58,60,61],:);
% anticancer (3, 4, 8, 10, 11, 15, 17, 18, 22, 24, 27, 29, 31, 33, 35, 37, 39, 40, 42) oxtrset7
% antimicrobial (45, 48, 52, 55, 59, 62, 63) oxtrset7

% cross examination study 4
%D4=D3([2,4,9,10,12,13,16,17,20,21,22,26,28,29,33,38,39,41,43,45,46,47,48,49,52,54,56,57,58,59,60,61],:);
% anticancer (3, 6, 7, 8, 11, 14, 15, 18, 24, 25, 27, 31, 32, 34, 35, 36, 37, 40, 42) oxtrset7
% antimicrobial (44, 50, 51, 53, 55, 62, 63) oxtrset7


sel=3;
%D4(sel,:)=D3(7,:);

%D=D4';
D=D3(2:63,:)';
[m,n]=size(D);

%figure(55), plot(1:21, D1(43,:),'ro');

%D=abs(D1);

Lambda1=0;
RR=zeros(n,n);
Z=(D'*D);
CC=zeros(m,n); Rc=zeros(m,n); HH=zeros(m,n); e=1; Vec=0;
p=1;
while e > 0
    if p > 1 
        Z=RR;
    end
    % calculating the eigenvectors and eigenvalues
    RR=Z-Lambda1*Vec*Vec';
    [Q,Lambda]=eig(RR);
    Lambda1=max(max(Lambda));
    if p == 1
        maxii=Lambda1;
    end
        
    for y=1:n
        for t=1:y
            if Lambda1 == Lambda(t,y);
                x=y;
            end
        end
    end
    Vec=Q(:,x);
    % Terminating the loop if the eigenvalue is less than 0.01 of the maximum eigenvalue
    if Lambda1 > maxii*0.14
        HH=D*Q;
        Rc(:,p)=HH(:,x);
        CC(p,:)=Vec';
    else
        e=0;
    end
    p=p+1;
end

%p=3; % for 2 important eigenvalues
p=4; % for 3 important eigenvalues
%p=5; % for 4 important eigenvalues
%p=6; % for 5 important eigenvalues
%p=7; % for 6 important eigenvalues
%p=8; % for 7 important eigenvalues

% Reconstructing the R & C matrices
R=zeros(m,p-1);
C=zeros(p-1,n);
for aa=1:p-1
    R(:,aa)=Rc(:,aa);
    C(aa,:)=CC(aa,:);
end
R;
C;

%per2factor=100*(Lambda(n,n)+Lambda(n-1,n-1))/sum(sum(Lambda));
%per3factor=100*(Lambda(n,n)+Lambda(n-1,n-1)+Lambda(n-2,n-2))/sum(sum(Lambda));
%per4factor=100*(Lambda(n,n)+Lambda(n-1,n-1)+Lambda(n-2,n-2)+Lambda(n-3,n-3))/sum(sum(Lambda));

per1factor=100*(abs(Lambda(1,1)))/sum(sum(abs(Lambda)));
per2factor=100*(abs(Lambda(1,1))+abs(Lambda(2,2)))/sum(sum(abs(Lambda)));
per3factor=100*(abs(Lambda(1,1))+abs(Lambda(2,2))+abs(Lambda(3,3)))/sum(sum(abs(Lambda)));
per4factor=100*(abs(Lambda(1,1))+abs(Lambda(2,2))+abs(Lambda(3,3))+abs(Lambda(4,4)))/sum(sum(abs(Lambda)));

for i=1:n
    pervar(1,i)=abs(100*Lambda(i,i)/sum(sum(Lambda)));
end

% Target Transformation
R1=R;
C1=C;
for phii=0:0
    phi=phii*pi/180;
%T=[cos(phi) -sin(phi); sin(phi) cos(phi)];
T=[1.0*cos(phi) -1.0*sin(phi) 0.0; 1.0*sin(phi) 1.0*cos(phi) 0.0; 0.0 0.0 1.0];
%T=eye(7);

R=R1*T;
C=inv(T)*C1;

figure(1)
plot(D(:,:),'-o'); axis([0 m+1 -3.5 3.5]);
title('Standardized Original Data (D)','FontSize', 18')
xlabel('QM Descriptors','FontSize', 18')
ylabel('Intensity','FontSize', 18')
set(gcf,'color','w');
set(gca,'FontSize',16);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontName','Arial','fontsize',8)

figure(2)
clf
plot(R*C(:,:),'-o'); axis([0 m+1 -3.5 3.5]);
title('Reconstructed Data with 3 PCs','FontSize', 18')
xlabel('QM Descriptors','FontSize', 18')
ylabel('Intensity','FontSize', 18')
set(gcf,'color','w');
set(gca,'FontSize',16);
%a = get(gca,'XTickLabel');
%set(gca,'XTickLabel',a,'FontName','Arial','fontsize',8)

figure(3)
plot(1:m, R(:,1),'ro-',1:m, R(:,2),'bx-', 1:m, R(:,3),'c+-');
title('Scores Matrix (R)','FontSize', 18')
xlabel('(Red: PC1, Blue: PC2, Cyan: PC3)','FontSize', 18')
ylabel('Intensity','FontSize', 18')
set(gcf,'color','w');
set(gca,'FontSize',16);
grid on
axis([0 16 -8.0 8.0]);

figure(4)
plot(1:n, C(1,:),'ro',1:n, C(2,:),'bx', 1:n, C(3,:),'c+');  %axis([-0.8 0.8 0 3]);
title('Loadings Matrix (C)','FontSize', 18')
xlabel('(Red: PC1, Blue: PC2, Cyan: PC3)','FontSize', 18')
ylabel('Intensity','FontSize', 18'); grid;
set(gcf,'color','w');
set(gca,'FontSize',16);
axis([0 n+1 -0.4 0.4]);

figure(5)
%plot(OxTrSet(2:29,1),C(1,:),'ro-',OxTrSet(2:29,1),C(2,:),'b+-',OxTrSet(2:29,1),C(3,:),'go-',OxTrSet(2:29,1),C(1,:)+C(2,:)+C(3,:),'yo-')
plot(1:n, abs(C(1,:)),'ro',1:n, abs(C(2,:)),'bx', 1:n, abs(C(1,:))+abs(C(2,:)),'go-');  %axis([-0.8 0.8 0 3]);
title('Loadings vs. actual')
xlabel('C1')
ylabel('Actual Loadings');
set(gcf,'color','w');

figure(6)
%plot(OxTrSet(2:29,1),C(1,:),'ro-',OxTrSet(2:29,1),C(2,:),'b+-',OxTrSet(2:29,1),C(3,:),'go-',OxTrSet(2:29,1),C(1,:)+C(2,:)+C(3,:),'yo-')
plot(1:n, (C(1,:)),'ro',1:n, (C(2,:)),'bx', 1:n, (C(3,:)),'c+', 1:n, (C(1,:)+C(2,:)+C(3,:)),'go-');  %axis([-0.8 0.8 0 3]);
title('Loadings vs. actual','FontSize', 16')
xlabel('C1','FontSize', 16')
ylabel('Actual Loadings','FontSize', 16'); grid;
set(gcf,'color','w');

%figure(7)
%plot(OxTrSet(2:29,1),C(1,:),'ro-',OxTrSet(2:29,1),C(2,:),'b+-',OxTrSet(2:29,1),C(3,:),'go-',OxTrSet(2:29,1),C(1,:)+C(2,:)+C(3,:),'yo-')
%plot(1:28, abs(C(1,:)),'ro',1:28, abs(C(2,:)),'bx', 1:28, abs(C(2,:))./abs(C(1,:)),'go-');  axis([0 30 0 1]);
%title('Loadings vs. actual')
%xlabel('C1')
%ylabel('Actual Loadings');
%set(gcf,'color','w');

figure(8)
%plot(OxTrSet(2:29,1),C(1,:),'ro-',OxTrSet(2:29,1),C(2,:),'b+-',OxTrSet(2:29,1),C(3,:),'go-',OxTrSet(2:29,1),C(1,:)+C(2,:)+C(3,:),'yo-')
plot(C(1,1:18), C(2,1:18),'ro',C(1,19:32), C(2,19:32),'b+');  %axis([-0.8 0.8 0 3]);
%plot(C(1,1:42), C(2,1:42),'ro',C(1,43:62), C(2,43:62),'b+');  %axis([-0.8 0.8 0 3]);
title('Red (Anticancer), Blue (Antimicrobial)','FontSize', 18')
xlabel('PC1 Loadings','FontSize', 18')
ylabel('PC2 Loadings','FontSize', 18'); grid;
set(gcf,'color','w');
set(gca,'FontSize',16);

figure(9)
%plot(OxTrSet(2:29,1),C(1,:),'ro-',OxTrSet(2:29,1),C(2,:),'b+-',OxTrSet(2:29,1),C(3,:),'go-',OxTrSet(2:29,1),C(1,:)+C(2,:)+C(3,:),'yo-')
plot(R(:,1), R(:,2),'ro');  %axis([-0.8 0.8 0 3]);
title('Loadings vs. actual')
xlabel('C1')
ylabel('Actual Loadings'); grid;
set(gcf,'color','w');
set(gca,'FontSize',14);

figure(10)
bar(1:n, pervar); axis([0.5 10.5 0 70]);
xlabel('Principal Components', 'FontSize', 18)
ylabel('Percentage of Explained Variances', 'FontSize',18);
set(gcf,'color','w');
set(gca,'FontSize',16);

%figure(11)
% plotted for manuscript 7 -> 3 study in cross examination set 1
%plot(1:n, C(1,:),'ro-',1:n, C(2,:),'bx-', 1:n, C(3,:),'c+-',1:n, CMatrix1814_2(1,:),'ro:',1:n, CMatrix1814_2(2,:),'bx:', 1:n, CMatrix1814_2(3,:),'c+:');  %axis([-0.8 0.8 0 3]);
%title('Loadings Matrix (C)','FontSize', 18')
%xlabel('(Red: PC1, Blue: PC2, Cyan: PC3)','FontSize', 18')
%ylabel('Intensity','FontSize', 18'); %grid;
%set(gcf,'color','w');
%set(gca,'FontSize',16);
%axis([0 33 -0.45 0.45]);

%figure(12)
% plotted for manuscript 7 -> 3 study in cross examination set 1
%plot(1:m, R(:,1),'ro-',1:m, R(:,2),'bx-', 1:m, R(:,3),'c+-',1:m, RMatrix1814_2(:,1),'ro:',1:m, RMatrix1814_2(:,2),'bx:', 1:m, RMatrix1814_2(:,3),'c+:');  %axis([-0.8 0.8 0 3]);
%title('Scores Matrix (R)','FontSize', 18')
%xlabel('(Red: PC1, Blue: PC2, Cyan: PC3)','FontSize', 18')
%ylabel('Intensity','FontSize', 18'); %grid;
%set(gcf,'color','w');
%set(gca,'FontSize',16);
%axis([0 33 -0.45 0.45]);

F=C';

if abs(F(1,1))-abs(F(2,1)) < 0.002
    a=1;
    b=2;
else
    a=2;
    b=1;
end

phii;
pause(0.5)
end

%save RMatrix1814_2.txt R -ascii -tabs -double;

sim=2;
if sign(C(1,sim))==sign(CMatrix_62(1,sim))
    C2(1,:)=C(1,:);
else
    C2(1,:)=-C(1,:);
end

if sign(C(2,sim))==sign(CMatrix_62(2,sim))
    C2(2,:)=C(2,:);
else
    C2(2,:)=-C(2,:);
end

if sign(C(3,sim))==sign(CMatrix_62(3,sim))
    C2(3,:)=C(3,:);
else
    C2(3,:)=-C(3,:);
end

for i=1:n
    CDiff(1:3,i)=CMatrix_62(1:3,i)-C2(1:3,sel);
    %CDiff(1:3,i)=abs(CMatrix1814(1:3,i))-abs(C(1:3,sel));
    %CDiff(1:3,i)=abs(CMatrixrot(1:3,i))-abs(C(1:3,sel));
end

CDiff2(1,:)=CDiff(1,:)*per1factor/100;
CDiff2(2,:)=CDiff(2,:)*(per2factor-per1factor)/100;
CDiff2(3,:)=CDiff(3,:)*(per3factor-per2factor)/100;
%CDiff2(4,:)=CDiff(4,:)*(per4factor-per3factor)/100;

AbsCDiff=abs(CDiff2);
AbsCDiff=sum(abs(CDiff2)); minAbsCDiff=min(AbsCDiff);
%AbsCDiff=abs(sum(CDiff)); minAbsCDiff=min(AbsCDiff);

AvgAbsCDiffCancer=sum(AbsCDiff(:,1:18))/18;
AvgAbsCDiffMicrobial=sum(AbsCDiff(:,19:32))/14;

if AvgAbsCDiffCancer < AvgAbsCDiffMicrobial
    disp('Anticancer based on average.');
else
    disp('Antimicrobial based on average.');
end

for jj=1:n
    if AbsCDiff(:,jj) == minAbsCDiff && jj < n-13
        disp('This compound is anticancer.');
        X = ['Similar to ',num2str(jj), ' with a minimum value of ', num2str(minAbsCDiff)];
        disp(X)
    elseif AbsCDiff(:,jj) == minAbsCDiff && jj >= n-13
        disp('This compound is antimicrobial.');
        X = ['Similar to ',num2str(jj), ' with a minimum value of ', num2str(minAbsCDiff)];
        disp(X)
    end
end

