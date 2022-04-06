function  Plotting_solar_stochastic(RSr,RSf,InitBody,v,x_road,roadall,lr,lf,theta,y2,y1f,y1r,time,timeall,Front_exct)

%Results magnification
% theta=5*theta;
% y2=5*y2;
% InitBody=5*InitBody;
y2f=y2+theta*lf;
y2r=y2-theta*lr;

%Road calculations
X_road=roadall-x_road+lf;
Y_road=Front_exct-0.5;

I = imread('solar.png'); %Read image 
BW = imbinarize(I); %PNG to binar
B = bwboundaries(BW,'noholes'); %Binar to boundaries
for k = 1:length(B)
    boundary = B{k};
end

R = [cos(-theta) -sin(-theta); sin(-theta) cos(-theta)];


X=-boundary(:,2)+800; 
X=X*0.00325;
Y=-boundary(:,1)+650;
Y=Y*0.00325;
BodyCord=[X,Y];
for i=1:length(X);
   BodyCord(i,:)=BodyCord(i,:)*R; 
end
X=BodyCord(:,1)-1.62+lf;
Y=BodyCord(:,2)+y2-InitBody;

figure (1); clf;
text(2.2,1.2,[num2str(time,'%2.1f') ' sec']); %Current time display
hold on
fill(X,Y,'w');              %Plotting Solar Body
hold on
scatter(0,y2,'filled','k');   %Plotting body GC
hold on
%Plotting wheels
r=0.3;
th = 0:pi/50:2*pi;
%Rear wheel
Rxunit = r * cos(th) - lr+0.01;
Ryunit = r * sin(th) - 0.2 + y1r;
plot(Rxunit, Ryunit,'r');
daspect([1 1 1])
hold on;
%Rear spring
Rdiff=-(y2r+0.2625-(-0.2+y1r));
line([-lr -lr-r/4 -lr+r/4 -lr-r/4 -lr+r/4 -lr-r/4 -lr+r/4 -lr-r/4 -lr+r/4 -lr-r/4 -lr],...
    [-0.2+y1r -0.2+y1r-Rdiff*0.5/10 -0.2+y1r-Rdiff*2/10 -0.2+y1r-Rdiff*3/10 ...
    -0.2+y1r-Rdiff*4/10 -0.2+y1r-Rdiff*5/10 -0.2+y1r-Rdiff*6/10 -0.2+y1r-Rdiff*7/10 ...
    -0.2+y1r-Rdiff*8/10 -0.2+y1r-Rdiff*9.5/10 y2r+0.2625],'Color','r');
%Front wheel
Fxunit = r * cos(th) + lf;
Fyunit = r * sin(th) - 0.2 + y1f;
plot(Fxunit, Fyunit,'b');
hold on;
%line([lf Fxunit],[-0.2+y1f Fyunit]);
%hold on;
%Front spring
Fdiff=-(y2f+0.2625-(-0.2+y1f));
line([lf lf-r/4 lf+r/4 lf-r/4 lf+r/4 lf-r/4 lf+r/4 lf-r/4 lf+r/4 lf-r/4 lf],...
    [-0.2+y1f -0.2+y1f-Fdiff*0.5/10 -0.2+y1f-Fdiff*2/10 -0.2+y1f-Fdiff*3/10 ...
    -0.2+y1f-Fdiff*4/10 -0.2+y1f-Fdiff*5/10 -0.2+y1f-Fdiff*6/10 -0.2+y1f-Fdiff*7/10 ...
    -0.2+y1f-Fdiff*8/10 -0.2+y1f-Fdiff*9.5/10 y2f+0.2625],'Color','b');
%Plotting road 
plot(X_road,Y_road,'k');
hold on;

%%line([lf(1) lf(1)-r/2 lf(1)+r/2 lf(1)],...
    %%[-0.2+y1f(1) -0.2+y1f(1)-RSf(1)*1/3 -0.2+y1f(1)-RSf(1)*2/3 y2f(1)]);

%axis([ min(X)-1 max(X)+1 min(Y)-1 max(Y)+1])
axis([ -4 3 -1 2])
daspect([1 1 1])
end