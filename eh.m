tStart=tic;
%Metricrelateddataandeqs.
G=1;%gravitationalconstant
M=1;%massofblackhole
a=0.6;%Kerr(spin)parameter
R=2∗G∗M;
radiuscelestialsphere=80;
aDiskMin=2∗R;
aDiskMax=5∗R;
%Plotblackspherefortheblackholeif’plottype
%andholdonthetheplot
ifplottype==1
[X,Y,Z]=sphere;
colormap([0,0,0])
surf(R∗X,R∗Y,R∗Z)
axisequal
holdon
end
=1’ischosen
%Somephysicalformulas(usedtolightennotation)
Sigma=@(r,theta)rˆ2+aˆ2∗cos(theta)ˆ2;
drSigma=@(r)2∗r;
dthetaSigma=@(theta)−2∗aˆ2∗cos(theta)∗sin(theta);
Delta=@(r)rˆ2−R∗r+aˆ2;
drDelta=@(r)2∗r−R;
%Physicaldimensionsofobservationalwindow
windowheight=0.00001;
windowwidth=(resolutionwidth/resolutionheight)∗windowheight;
distancefromwindow=−2∗0.000007;
%Initializecoordinatesmatricesforimagemapping
%coordsnoaDisk(i,j,1):theta
%coordsnoaDisk(i,j,2):phi
%coordsnoaDisk(i,j,3):1ifstoppedatR
%
0ifstoppedatradiuscelestialsphere
coordsnoaDisk=zeros(resolutionheight,resolutionwidth,3);
%coordsaDisk(i,j,1):radius
%coordsaDisk(i,j,2):phi
%coordsaDisk(i,j,3):0ifdidnothitaccretiondisk
%
1ifpassthroughaccretiondisk
coordsaDisk=zeros(resolutionheight,resolutionwidth,3);
stepsize=0.1;%stepsizeforRunge−Kutta4
hbar=parforprogressbar(resolutionwidth,’Pleasewait...’);%createtheprogressbar(resolutionwidth);
forj=1:resolutionwidth
hbar.iterate(1);%updateprogressbyoneiteration
fori=1:resolutionheight
%
%
%Pixelslocationonobservationalwindow(weommitthevertical
%singularityatpi/2byintroducingajump:singularityhack)
%singularityhack=0.01;
ifj<resolutionwidth/2−singularityhack∗resolutionwidth
h=windowheight/2−(i−1)∗windowheight/(resolutionheight−1);
w=−windowwidth/2+(j−1)∗windowwidth/(resolutionwidth−1);
else
h=windowheight/2−(i−1)∗windowheight/(resolutionheight−1);
w=−windowwidth/2+(singularityhack∗windowwidth)+(j−1)∗windowwidth/(resolutionwidth−1);
end
h=windowheight/2−(i−1)∗windowheight/(resolutionheight−1);
w=−windowwidth/2+(j−1)∗windowwidth/(resolutionwidth−1);
%Initializinginitialconditions
r=70;
theta=pi/2−pi/46;%offsettheblackholetoseewithofaDisk
phi=0;
tdot=1;
phidot=(csc(theta)∗w)/sqrt((aˆ2+rˆ2)∗(distancefromwindowˆ2+wˆ2+hˆ2));
pr=2∗Sigma(r,theta)∗(h∗(aˆ2+rˆ2)∗cos(theta)+r∗sqrt(aˆ2+rˆ2)∗sin(theta)∗distancefromwindow
/(sqrt(distancefromwindowˆ2+hˆ2+wˆ2)∗(aˆ2+2∗rˆ2+aˆ2∗cos(2∗theta))∗Delta(r));
ptheta=2∗Sigma(r,theta)∗(−h∗r∗sin(theta)+sqrt(aˆ2+rˆ2)∗cos(theta)∗distancefromwindow)...
/(sqrt(distancefromwindowˆ2+hˆ2+wˆ2)∗(aˆ2+2∗rˆ2+aˆ2∗cos(2∗theta)));
)...
%Conservedquantities(instantiatedwithI.C.)
E=(1−R/r)∗tdot+(R∗a∗phidot)/r;
L=−(R∗a)/r∗tdot+(rˆ2+aˆ2+(R∗aˆ2)/r)∗phidot;
%Geodesicequations(systemoffirstorderODEs)
%input:x=[r;theta;phi;pr;ptheta]
%ouput:dx=[dr;dtheta;dphi;dpr;dptheta]”xdot”
f=@(lambda,x)
[(x(4)∗Delta(x(1)))/Sigma(x(1),x(2));
x(5)/Sigma(x(1),x(2));
(a∗(−a∗L+x(1)∗R∗E)+L∗csc(x(2))ˆ2∗Delta(x(1)))/(Delta(x(1))∗Sigma(x(1),x(2)));
−(1/(2∗Delta(x(1))ˆ2∗Sigma(x(1),x(2))ˆ2))∗(Sigma(x(1),x(2))∗(−E∗Delta(x(1))∗...
(a∗R∗(−2∗L+a∗E∗sin(x(2))ˆ2)+2∗x(1)∗E∗Sigma(x(1),x(2)))...
+(a∗(a∗Lˆ2−2∗L∗x(1)∗R∗E+a∗x(1)∗R∗Eˆ2∗sin(x(2))ˆ2)+x(4)ˆ2∗...
Delta(x(1))ˆ2+(aˆ2+x(1)ˆ2)∗Eˆ2∗Sigma(x(1),x(2)))∗drDelta(x(1)))...
+Delta(x(1))∗(a∗(L∗(a∗L−2∗x(1)∗R∗E)+a∗x(1)∗R∗Eˆ2∗sin(x(2))ˆ2)...
−Delta(x(1))∗(x(5)ˆ2+Lˆ2∗csc(x(2))ˆ2+x(4)ˆ2∗Delta(x(1))))∗drSigma(x(1)));
−(1/(2∗Delta(x(1))∗Sigma(x(1),x(2))ˆ2))∗(−2∗sin(x(2))∗(aˆ2∗x(1)∗R∗Eˆ2∗cos(x(2))...
+Lˆ2∗cot(x(2))∗csc(x(2))ˆ3∗Delta(x(1)))∗Sigma(x(1),x(2))...
+(a∗(L∗(a∗L−2∗x(1)∗R∗E)+a∗x(1)∗R∗Eˆ2∗sin(x(2))ˆ2)−Delta(x(1))...
∗(x(5)ˆ2+Lˆ2∗csc(x(2))ˆ2+x(4)ˆ2∗Delta(x(1))))∗dthetaSigma(x(2)))];
%SolvingforthegeodesicsusingRunga−Kutta4
%RK4parameters
x0=[r;theta;phi;pr;ptheta];
curve=x0’;
switchplottype
case1
k=1;
%Curvesarecomputedwithinthecelestialspherewheretheyexist
while(R<curve(k,1))&&(curve(k,1)<radiuscelestialsphere)&&(k<20000)
%Cleancoordintatesvalues
curve(k,2)=mod(curve(k,2),2∗pi);
curve(k,3)=mod(curve(k,3),2∗pi);
ifcurve(k,2)>pi
curve(k,2)=2∗pi−curve(k,2);
curve(k,3)=mod(pi+curve(k,3),2∗pi);
end
k=k+1;
%UseRunge−Kutta4steptoevolvegeodesiccurve
%Someadaptativityisobtainedbymultiplyingstepsize
%byDelta,whichscaleinverslytothecoordinates
%singularityneartheeventhorizon
curve(k,:)=rk4step(f,0,curve(k−1,:)’,min([stepsize∗Delta(curve(k−1,1));stepsize]))’;
end
%Transformtoeuclideancoordinatesandplot3
[n,m]=size(curve);
A=a∗ones(n,1);
Boyer2Cart=@(r,theta,phi)[sqrt(r.ˆ2+A.ˆ2).∗sin(theta).∗cos(phi),...
sqrt(r.ˆ2+A.ˆ2).∗sin(theta).∗sin(phi),...
r.∗cos(theta)];
cart=Boyer2Cart(curve(:,1),curve(:,2),curve(:,3));
plot3(cart(:,1),cart(:,2),cart(:,3));
case0
k=1;%keepingcountofsteps
passedthroughaDisk=0;
%Curvesstepsarecomputedwithinthecelestialspherewheretheyexist
while(1.2∗R<curve(1))&&(curve(1)<radiuscelestialsphere)&&(k<20000)
%UseRunge−Kuttatotakeatemporarystepfoward
tempcurvestep=rk4step(f,0,curve’,min([stepsize∗Delta(curve(1));stepsize]))’;
ifpassedthroughaDisk==0
%Checkifthetemporarystepgothroughtheaccretiondisk
if(tempcurvestep(2)−pi/2)∗(curve(2)−pi/2)<0
if(aDiskMin<tempcurvestep(1))&&(tempcurvestep(1)<aDiskMax)
1];
coordsaDisk(i,j,:)=[tempcurvestep(1);tempcurvestep(3);
passedthroughaDisk=1;
end
end
end
%Updatecurvewiththetemporarystep
curve=tempcurvestep;
%Cleancoordinatesvalues
curve(2)=mod(curve(2),2∗pi);
curve(3)=mod(curve(3),2∗pi);
ifcurve(2)>pi
curve(2)=2∗pi−curve(2);
curve(3)=mod(pi+curve(3),2∗pi);
end
k=k+1;
end
%Savecoordinatesaccordingtocases
if(1.2∗R<curve(1))
coordsnoaDisk(i,j,:)=[curve(2);curve(3);0];
else
coordsnoaDisk(i,j,:)=[curve(2);curve(3);1];
end
end
end
end
close(hbar);
ifplottype==1
holdoff
end
%Savecoordinatesasmaltabfile
saveimgMapnoBHcoordsnoaDiskcoordsaDiskaDiskMinaDiskMax;
%%Generatingfinalimagesfromdifferent
celestialscenes
imagetypewithdisk=1;
booleanaDisk=1;
load(’imgMapnoBH’);
[resolutionheight,resolutionwidth,˜]=size(coordsnoaDisk);
accretiondiskscene=imread(’adiskskewed.png’);
[accretionsceneheightres,accretionscenewidthres,˜]=size(accretiondiskscene);
k=1;
switchk
case1
celestialscene=imread(’GriddedGal.jpg’);
case2
celestialscene=imread(’InterstellarWormholeFig6b.jpg’);
case3
celestialscene=imread(’InterstellarWormholeFig10.jpg’);
end
%Resolutiondatarelatedtoscene
[celestialsceneheightres,celestialscenewidthres,˜]=size(celestialscene);
IMG=zeros(resolutionheight,resolutionwidth,3,’uint8’);
%Generateimage
forj=1:resolutionwidth
fori=1:resolutionheight
ifimagetypewithdisk==1
ifbooleanaDisk==1
ifcoordsaDisk(i,j,3)==1
%Cleancoordinatesvalues
coordsaDisk(i,j,2)=mod(coordsaDisk(i,j,2),2∗pi);
ifcoordsaDisk(i,j,2)>pi
coordsaDisk(i,j,2)=2∗pi−coordsaDisk(i,j,2);
end
%coordsaDisk(i,j,2)
%round(coordsaDisk(i,j,2)∗((accretionscenewidthres−1)/(2∗pi))+1)
ifcoordsnoaDisk(i,j,3)==1
IMG(i,j,:)=accretiondiskscene(round((coordsaDisk(i,j,1)−aDiskMin)∗((accretionsceneheightres−1)/(aDiskMax−aDiskMin))+1),...
round(coordsaDisk(i,j,2)∗((accretionscenewidthres−1)/(2∗pi))+1),:);
else
IMG(i,j,:)=max(accretiondiskscene(round((coordsaDisk(i,j,1)−aDiskMin)∗((accretionsceneheightres−1)/(aDiskMax−aDiskMin))+1),...
round(coordsaDisk(i,j,2)∗((accretionscenewidthres−1)/(2∗pi))+1),:),...
celestialscene(round(coordsnoaDisk(i,j,1)∗((celestialsceneheightres−1)/pi)+1),...
round(coordsnoaDisk(i,j,2)∗((celestialscenewidthres−1)/(2∗pi))+1),:));
end
booleanaDisk=0;
end
end
end
ifbooleanaDisk==1
ifcoordsnoaDisk(i,j,3)==1
IMG(i,j,1)=0;
IMG(i,j,2)=0;
IMG(i,j,3)=0;
else
IMG(i,j,:)=celestialscene(round(coordsnoaDisk(i,j,1)∗((celestialsceneheightres−1)/pi)+1),...
round(coordsnoaDisk(i,j,2)∗((celestialscenewidthres−1)/(2∗pi))+1),:);
end
end
booleanaDisk=1;
end
end
IMG=[IMG(:,1:(resolutionwidth/2−3),:)IMG(:,(resolutionwidth/2−3),:)max(IMG(:,(resolutionwidth/2−3),:),IMG(:,resolutionwidth/2+5,:))...
IMG(:,resolutionwidth/2+5,:)IMG(:,(resolutionwidth/2+5):resolutionwidth,:)];
switchk
case1
imwrite(IMG,’grid’,’jpg’)
case2
imwrite(IMG,’test2InterstellarWormholeFig6BHdisk.jpg’,’jpg’)
case3
imwrite(IMG,’test2InterstellarWormholeFig10BHdisk.jpg’,’jpg’)
end
%Finalcomputationstime
tEnd=toc(tStart);
fprintf(’%dminutesand%fseconds\n’,floor(tEnd/60),rem(tEnd,60));