clear all;
clc;
format long;
E=2e11; %Young's modulus
G=7.692307692307692e10; %shear modulus
k=5/6; %section shape factor (rectangular cross-section)
I=0.05*0.02^3/12; %Second moment of area
Ip=5.34166667e-7; %polar second moment of area
A=0.05*0.02; %Cross-sectional area
Rho=7850; %density
b=0.05; %Cross-sectional width
h=0.02; %Cross-sectional height
nodes=8;
conn=[1 2; 2 3; 3 4; 4 5; 5 6; 6 7; 7 8]; %matrix of the two nodes connecting each member, the crack element is regarded as a zero-length member 
angle=[pi/2 pi/2 pi/2 0 0 0 -pi/2];
sup1=[1 2 3]; %left boundary condition
sup2=[22 23 24]; %right boundary condition
l=1; %length of the beam member
Xi1=0.2; %Crack depth to sectional height ratio
Xi2=0.3;
Mu=0.3; %Poisson's ratio
F1(1,1)=exp(1)^(1/(1-Xi1))*(-0.326584e-5*Xi1+1.455190*Xi1^2-0.984690*Xi1^3+4.895396*Xi1^4-6.501832*Xi1^5+12.792091*Xi1^6-26.723556*Xi1^7+35.073593*Xi1^8-34.954632*Xi1^9+9.054062*Xi1^10);
F1(2,2)=exp(1)^(1/(1-Xi1))*(-0.326018e-6*Xi1+1.454954*Xi1^2-1.455784*Xi1^3-0.421981*Xi1^4-0.279522*Xi1^5+0.455399*Xi1^6-2.432830*Xi1^7+5.427219*Xi1^8-6.643057*Xi1^9+4.466758*Xi1^10);
F1(3,3)=exp(1)^(1/(1-Xi1))*(-0.219628e-4*Xi1+52.37903*Xi1^2-130.2483*Xi1^3+308.442769*Xi1^4-602.445544*Xi1^5+939.044538*Xi1^6-1310.95029*Xi1^7+1406.52368*Xi1^8-1067.4998*Xi1^9+391.536356*Xi1^10);
F2(1,1)=exp(1)^(1/(1-Xi2))*(-0.326584e-5*Xi2+1.455190*Xi2^2-0.984690*Xi2^3+4.895396*Xi2^4-6.501832*Xi2^5+12.792091*Xi2^6-26.723556*Xi2^7+35.073593*Xi2^8-34.954632*Xi2^9+9.054062*Xi2^10);
F2(2,2)=exp(1)^(1/(1-Xi2))*(-0.326018e-6*Xi2+1.454954*Xi2^2-1.455784*Xi2^3-0.421981*Xi2^4-0.279522*Xi2^5+0.455399*Xi2^6-2.432830*Xi2^7+5.427219*Xi2^8-6.643057*Xi2^9+4.466758*Xi2^10);
F2(3,3)=exp(1)^(1/(1-Xi2))*(-0.219628e-4*Xi2+52.37903*Xi2^2-130.2483*Xi2^3+308.442769*Xi2^4-602.445544*Xi2^5+939.044538*Xi2^6-1310.95029*Xi2^7+1406.52368*Xi2^8-1067.4998*Xi2^9+391.536356*Xi2^10);
Lcr11=(1-Mu^2)/(E*b)*F1(1,1);
Lcr22=(1-Mu^2)/(E*b)*F1(2,2);
Lcr33=(1-Mu^2)/(E*b*h^2)*F1(3,3);
Mcr11=(1-Mu^2)/(E*b)*F2(1,1);
Mcr22=(1-Mu^2)/(E*b)*F2(2,2);
Mcr33=(1-Mu^2)/(E*b*h^2)*F2(3,3);
%%%%%% data for mode shape calculation %%%%%%
nodesdiv=303; %The number of nodes dividing the cracked frame
Ldiv([1:50 52:121 123:302])=0.01; %the length of beam elements %the crack is regarded as zero-length elements
coorx(1:102)=0;
coorx(102:122)=0:0.01:0.2;
coorx(123:203)=0.2:0.01:1;
coorx(203:303)=1;
coory(1:51)=0:0.01:0.5;
coory(52:102)=0.5:0.01:1;
coory(102:203)=1;
coory(203:303)=1:-0.01:0;
% X and Y coordinates of the small elements' nodes, required for mode shape plotting
numelem=302;
for i=1:numelem
	connect(i,:)=[i i+1];
end
angles(1:101)=pi/2;
angles(102:202)=0;
angles(203:302)=-pi/2;
sup11=[1 2 3];
sup22=[907 908 909];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reqmodes=[1 2 3 4 5 6 7 8 9 10]; %Required natural frequencies, must be arranged in ascending order of modes
fr=length(reqmodes); %Total number of required natural frequencies
L(1)=0.5;
L(2)=0;
L(3)=0.5;
L(4)=0.2;
L(5)=0;
L(6)=0.8;
L(7)=1;
w=0.0001; %Trial value for 'w'
J0=frame_WWalgorithm(w,L,nodes,conn,angle,sup1,sup2,E,G,k,I,Ip,A,Rho,Mu,Lcr11,Lcr22,Lcr33,Mcr11,Mcr22,Mcr33);
for r=reqmodes %Calculating the required natural frequencies, the orders of which are defined by 'reqmodes'
    while r>J0 %Establishing the lower and upper limits around the required natural frequency
        wl=w; w=2*w; wu=w;
		J0=frame_WWalgorithm(w,L,nodes,conn,angle,sup1,sup2,E,G,k,I,Ip,A,Rho,Mu,Lcr11,Lcr22,Lcr33,Mcr11,Mcr22,Mcr33);
	end
	while (wu-wl)>0.0000001 %Required accuracy
	    w=(wu+wl)/2; %Bisection
		J0=frame_WWalgorithm(w,L,nodes,conn,angle,sup1,sup2,E,G,k,I,Ip,A,Rho,Mu,Lcr11,Lcr22,Lcr33,Mcr11,Mcr22,Mcr33);
		if r>J0
		    wl=w;
		else
		    wu=w;
		end
	end
	wr(r)=w;
	modeshape(r,:)=frame_WWmodeshape(w,connect,angles,sup11,sup22,nodesdiv,Ldiv,E,G,k,I,Ip,A,Rho,Mu,Lcr11,Lcr22,Lcr33,Mcr11,Mcr22,Mcr33);
end
fHz=wr./(2*pi); %Natural frequencies in Hz
for mode=reqmodes; %Mode shape calculations and plotting
    dx([1 303])=0;
	dy([1 303])=0; %clamped nodes
	dx(2:302)=0.04.*modeshape(mode,1:3:901);
	dy(2:302)=0.04.*modeshape(mode,2:3:902);
	figure
	plot(coorx(1:303)+dx(1:303), coory(1:303)+dy(1:303),'black','linewidth',3)
end