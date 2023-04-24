function[Modeshapes]=frame_WWmodeshape(w,connect,angles,sup11,sup22,nodesdiv,Ldiv,E,G,k,I,Ip,A,Rho,Mu,Lcr11,Lcr22,Lcr33,Mcr11,Mcr22,Mcr33)
%'w' is the trial frequency in rad/s.
DFGstiff=zeros(3*nodesdiv,3*nodesdiv); %Global dynamic stiffness matrix of the cracked frame
DFstiff=zeros(6,6); %dynamic stiffness matrix of each element, including the crack element.
Lcrackstiff=[1/Lcr11 0 0 -1/Lcr11 0 0; 0 1/Lcr22 0 0 -1/Lcr22 0; 0 0 1/Lcr33 0 0 -1/Lcr33; -1/Lcr11 0 0 1/Lcr11 0 0; 0 -1/Lcr22 0 0 1/Lcr22 0; 0 0 -1/Lcr33 0 0 1/Lcr33]; 
Mcrackstiff=[1/Mcr11 0 0 -1/Mcr11 0 0; 0 1/Mcr22 0 0 -1/Mcr22 0; 0 0 1/Mcr33 0 0 -1/Mcr33; -1/Mcr11 0 0 1/Mcr11 0 0; 0 -1/Mcr22 0 0 1/Mcr22 0; 0 0 -1/Mcr33 0 0 1/Mcr33];
for i=1:50
	Alpha2(i)=Rho*A*w^2*Ldiv(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*Ldiv(i)^4/(E*I);
	rr(i)=I/(A*Ldiv(i)^2);
	ss(i)=E*I/(k*A*G*Ldiv(i)^2);
	a1(i)=E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	S(i)=sin(Phi(i));
	C(i)=cos(Phi(i));
	W1(i)=E*I/Ldiv(i);
	W2(i)=E*I/Ldiv(i)^2;
	W3(i)=E*I/Ldiv(i)^3;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	if bb(i)*rr(i)*ss(i)<1
		j=1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sinh(Lambda(i));
		C_(i)=cosh(Lambda(i));	
	end
	if bb(i)*rr(i)*ss(i)>1
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sin(Lambda(i));
		C_(i)=cos(Lambda(i));	
	end
	Eta(i)=Z(i)/(j*Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-C(i)*C_(i))+(1-j*Eta(i)^2)*S(i)*S_(i));
	d1(i)=W3(i)*bb(i)*Tau(i)*(C(i)*S_(i)+Eta(i)*S(i)*C_(i))/(Lambda(i)*Phi(i));
	d2(i)=W2(i)*Z(i)*Tau(i)*((Phi(i)+j*Eta(i)*Lambda(i))*S(i)*S_(i)-(Lambda(i)-Eta(i)*Phi(i))*(1-C(i)*C_(i)))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=W1(i)*Tau(i)*(S(i)*C_(i)-j*Eta(i)*C(i)*S_(i));
	d4(i)=-W3(i)*bb(i)*Tau(i)*(S_(i)+Eta(i)*S(i))/(Lambda(i)*Phi(i));
	d5(i)=W2(i)*Z(i)*Tau(i)*(C_(i)-C(i));
	d6(i)=W1(i)*Tau(i)*(j*Eta(i)*S_(i)-S(i));
	DFstiff=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
	t=[cos(angles(i)) sin(angles(i)) 0; -sin(angles(i)) cos(angles(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=connect(i,1);
	n2=connect(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=51
	DFstiff=Lcrackstiff;
	t=[cos(angles(i)) sin(angles(i)) 0; -sin(angles(i)) cos(angles(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=connect(i,1);
	n2=connect(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=52:121
	Alpha2(i)=Rho*A*w^2*Ldiv(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*Ldiv(i)^4/(E*I);
	rr(i)=I/(A*Ldiv(i)^2);
	ss(i)=E*I/(k*A*G*Ldiv(i)^2);
	a1(i)=E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	S(i)=sin(Phi(i));
	C(i)=cos(Phi(i));
	W1(i)=E*I/Ldiv(i);
	W2(i)=E*I/Ldiv(i)^2;
	W3(i)=E*I/Ldiv(i)^3;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	if bb(i)*rr(i)*ss(i)<1
		j=1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sinh(Lambda(i));
		C_(i)=cosh(Lambda(i));	
	end
	if bb(i)*rr(i)*ss(i)>1
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sin(Lambda(i));
		C_(i)=cos(Lambda(i));	
	end
	Eta(i)=Z(i)/(j*Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-C(i)*C_(i))+(1-j*Eta(i)^2)*S(i)*S_(i));
	d1(i)=W3(i)*bb(i)*Tau(i)*(C(i)*S_(i)+Eta(i)*S(i)*C_(i))/(Lambda(i)*Phi(i));
	d2(i)=W2(i)*Z(i)*Tau(i)*((Phi(i)+j*Eta(i)*Lambda(i))*S(i)*S_(i)-(Lambda(i)-Eta(i)*Phi(i))*(1-C(i)*C_(i)))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=W1(i)*Tau(i)*(S(i)*C_(i)-j*Eta(i)*C(i)*S_(i));
	d4(i)=-W3(i)*bb(i)*Tau(i)*(S_(i)+Eta(i)*S(i))/(Lambda(i)*Phi(i));
	d5(i)=W2(i)*Z(i)*Tau(i)*(C_(i)-C(i));
	d6(i)=W1(i)*Tau(i)*(j*Eta(i)*S_(i)-S(i));
	DFstiff=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
	t=[cos(angles(i)) sin(angles(i)) 0; -sin(angles(i)) cos(angles(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=connect(i,1);
	n2=connect(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=122
	DFstiff=Mcrackstiff;
	t=[cos(angles(i)) sin(angles(i)) 0; -sin(angles(i)) cos(angles(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=connect(i,1);
	n2=connect(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
for i=123:302
	Alpha2(i)=Rho*A*w^2*Ldiv(i)^2/(E*A);
	Beta2(i)=Rho*Ip*Mu^2*w^2/(E*A);
	Gamma(i)=(Alpha2(i)/(1-Beta2(i)))^0.5;
	bb(i)=Rho*A*w^2*Ldiv(i)^4/(E*I);
	rr(i)=I/(A*Ldiv(i)^2);
	ss(i)=E*I/(k*A*G*Ldiv(i)^2);
	a1(i)=E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*cot(Gamma(i));
	a2(i)=-E*A/Ldiv(i)*Gamma(i)*(1-Beta2(i))*csc(Gamma(i));
	Phi(i)=(bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5)^0.5;
	S(i)=sin(Phi(i));
	C(i)=cos(Phi(i));
	W1(i)=E*I/Ldiv(i);
	W2(i)=E*I/Ldiv(i)^2;
	W3(i)=E*I/Ldiv(i)^3;
	Z(i)=Phi(i)-bb(i)*ss(i)/Phi(i);
	if bb(i)*rr(i)*ss(i)<1
		j=1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sinh(Lambda(i));
		C_(i)=cosh(Lambda(i));	
	end
	if bb(i)*rr(i)*ss(i)>1
		j=-1;
		Lambda(i)=(j*(-bb(i)*(rr(i)+ss(i))/2+bb(i)/2*((rr(i)+ss(i))^2+4/bb(i)*(1-bb(i)*rr(i)*ss(i)))^0.5))^0.5;
		S_(i)=sin(Lambda(i));
		C_(i)=cos(Lambda(i));	
	end
	Eta(i)=Z(i)/(j*Lambda(i)+bb(i)*ss(i)/Lambda(i));
	Tau(i)=(Lambda(i)+Eta(i)*Phi(i))/(2*Eta(i)*(1-C(i)*C_(i))+(1-j*Eta(i)^2)*S(i)*S_(i));
	d1(i)=W3(i)*bb(i)*Tau(i)*(C(i)*S_(i)+Eta(i)*S(i)*C_(i))/(Lambda(i)*Phi(i));
	d2(i)=W2(i)*Z(i)*Tau(i)*((Phi(i)+j*Eta(i)*Lambda(i))*S(i)*S_(i)-(Lambda(i)-Eta(i)*Phi(i))*(1-C(i)*C_(i)))/(Lambda(i)+Eta(i)*Phi(i));
	d3(i)=W1(i)*Tau(i)*(S(i)*C_(i)-j*Eta(i)*C(i)*S_(i));
	d4(i)=-W3(i)*bb(i)*Tau(i)*(S_(i)+Eta(i)*S(i))/(Lambda(i)*Phi(i));
	d5(i)=W2(i)*Z(i)*Tau(i)*(C_(i)-C(i));
	d6(i)=W1(i)*Tau(i)*(j*Eta(i)*S_(i)-S(i));
	DFstiff=[a1(i) 0 0 a2(i) 0 0; 0 d1(i) d2(i) 0 d4(i) d5(i); 0 d2(i) d3(i) 0 -d5(i) d6(i); a2(i) 0 0 a1(i) 0 0; 0 d4(i) -d5(i) 0 d1(i) -d2(i); 0 d5(i) d6(i) 0 -d2(i) d3(i)];
	t=[cos(angles(i)) sin(angles(i)) 0; -sin(angles(i)) cos(angles(i)) 0; 0 0 1];
	T=[t zeros(3,3); zeros(3,3) t];
	DFstiff=T'*DFstiff*T;
	n1=connect(i,1);
	n2=connect(i,2);
	DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)=DFGstiff(3*n1-2:3*n1,3*n1-2:3*n1)+DFstiff(1:3,1:3);
	DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)=DFGstiff(3*n1-2:3*n1,3*n2-2:3*n2)+DFstiff(1:3,4:6);
	DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)=DFGstiff(3*n2-2:3*n2,3*n1-2:3*n1)+DFstiff(4:6,1:3);
	DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)=DFGstiff(3*n2-2:3*n2,3*n2-2:3*n2)+DFstiff(4:6,4:6);
end
if isempty(sup11)&&isempty(sup22) %In the case of no supports or restraints.
    dis(1:(3*nodesdiv-1))=(inv(DFGstiff(1:(3*nodesdiv-1),1:(3*nodesdiv-1))))*DFGstiff(1:(3*nodesdiv-1),3*nodesdiv);
	dis(3*nodesdiv)=-1; %The highest numbered degree of freedom is assumed to have a unit displacement (or rotation)
	Modeshapes=dis;
else
    DFGstiff([sup22],:)=[];
    DFGstiff(:,[sup22])=[];
	DFGstiff([sup11],:)=[];
	DFGstiff(:,[sup11])=[];%Delete the rows and columns in the stiffness matrix, corresponding to the suppressed degrees of freedom.
    dis(1:(3*nodesdiv-length(sup22)-length(sup11)-1))=(inv(DFGstiff(1:(3*nodesdiv-length(sup22)-length(sup11)-1),1:(3*nodesdiv-length(sup22)-length(sup11)-1))))*DFGstiff(1:(3*nodesdiv-length(sup22)-length(sup11)-1),(3*nodesdiv-length(sup22)-length(sup11)));
	dis((3*nodesdiv-length(sup22)-length(sup11)))=-1;
	Modeshapes=dis;
end
end