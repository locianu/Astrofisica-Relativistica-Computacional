%RK-4 PARA PAPER RENAULD
%CONDI��ES INICIAIS

plot_type=1;

resolution_height=2; %resolu��es da tela, em pixels
resolution_width=2;

tStart=tic;

G=1;
M=1;  %massa do buraco negro
a=0.6;   %momento angular
R=2*G*M;
radius_celestial_sphere=80;
aDiskMin=2*R;
aDiskMax=5*R;

Sigma=@(r,theta)r^2+a^2*cos(theta)^2;
drSigma=@(r)2*r;
dthetaSigma=@(theta)-2*a^2*cos(theta)*sin(theta);
Delta=@(r)r^2-R*r+a^2;
drDelta=@(r)2*r-R;

%Plot black sphere for the blackhole if �plot_type=1� is chosen
%and hold on the plot
if plot_type==1
    [X,Y,Z]=sphere;
    colormap([0,0,0]);
    surf(R*X,R*Y,R*Z);
    axis equal;
    hold on;
end

%dimens�es da janela observacional
window_height=0.00001;
window_width=(resolution_width/resolution_height)*window_height;
distance_from_window=-1.4e-4;

coords_no_aDisk=zeros(resolution_height,resolution_width,3);
coords_aDisk=zeros(resolution_height,resolution_width,3);

stepsize=0.1; %passo do rk4






%hbar=parfor_progressbar(resolution_width,'Please wait ...'); %create the progress bar(resolution_width)%
fig=1;
for j=1:resolution_width
%    hbar.iterate(1);%update progress by one iteration
   for i=1:resolution_height
% i=4;
% j=4;
        h=window_height/2-(i-1)*window_height/(resolution_height-1);
        w=-window_width/2+(j-1)*window_width/(resolution_width-1);

        r=70;
        theta=pi/2-pi/46;%offset the blackhole to see with of aDisk
        phi=0;
        t_dot=1;

        phi_dot=(csc(theta)*w)/sqrt((a^2+r^2)*(distance_from_window^2+w^2+h^2));

        p_r=2*Sigma(r,theta)*(h*(a^2+r^2)*cos(theta)+r.*sqrt(a^2+r^2)*sin(theta)*distance_from_window)...
            /(sqrt(distance_from_window^2+h^2+w^2)*(a^2+2*r^2+a^2*cos(2*theta))*Delta(r));

        p_theta=2*Sigma(r,theta)*(-h*r*sin(theta)+sqrt(a^2+r^2)*cos(theta)*distance_from_window)...
            /(sqrt(distance_from_window^2+h^2+w^2)*(a^2+2*r^2+a^2*cos(2*theta)));

        E=(1-R/r)*t_dot+(R*a*phi_dot)/r;
        L=-(R*a)/r*t_dot+(r^2+a^2+(R*a^2)/r)*phi_dot;


        %GEOD�SICAS
        f1=@(r,theta,p_r,p_theta)[(p_r*Delta(r))/Sigma(r,theta)];

        f2=@(r,theta,p_r,p_theta)[(p_theta)/Sigma(r,theta)];

        f3=@(r,theta,p_r,p_theta)[-(1/(2*Delta(r)^2*Sigma(r,theta)^2))*(Sigma(r,theta)*...
            (-E*Delta(r)*(a*R*(-2*L+a*E*sin(theta)^2)+2*r*E*Sigma(r,theta))...
            +(a*(a*L^2-2*L*r*R*E+a*r*R*E^2*sin(theta)^2)+p_r^2*...
            Delta(r)^2+(a^2+r^2)*E^2*Sigma(r,theta))*drDelta(r))...
            +Delta(r)*(a*(L*(a*L-2*r*R*E)+a*r*R*E^2*sin(theta)^2)...
            -Delta(r)*(p_theta^2+L^2*csc(theta)^2+p_r^2*Delta(r)))*drSigma(r))];

        f4=@(r,theta,p_r,p_theta)[-(1/(2*Delta(r)*Sigma(r,theta)^2))*(-2*sin(theta)*(a^2*r*R*E^2*cos(theta)...
            +L^2*cot(theta)*csc(theta)^3*Delta(r))*Sigma(r,theta)...
            +(a*(L*(a*L-2*r*R*E)+a*r*R*E^2*sin(theta)^2)-Delta(r)...
            *(p_theta^2+L^2*csc(theta)^2+p_r^2*Delta(r)))*dthetaSigma(theta))];

        f5=@(r,theta,p_r,p_theta)(a*(-a*L+r*R*E)+L*csc(theta)^2*Delta(r))/(Delta(r)*Sigma(r,theta));

        x_0 = [r theta p_r p_theta phi];
        curve=x_0;
        %disp(x_0);
        %disp(curve);


        k=1
        Nk=20000

        while((R<r)&&(r<radius_celestial_sphere)&&(k<Nk))
            %Clean coordinates values
            curve(k,2)=mod(curve(k,2),2*pi);
            curve(k,5)=mod(curve(k,5),2*pi);
            if curve(k,2)>pi
                curve(k,2)=2*pi-curve(k,2);
                curve(k,5)=mod(pi+curve(k,5),2*pi);
            end
            theta = curve(k,2)

            %runge-kutta
            step=min([stepsize*Delta(r);stepsize])
            k1=step*f1(r,theta,p_r,p_theta);
            m1=step*f2(r,theta,p_r,p_theta);
            n1=step*f3(r,theta,p_r,p_theta);
            s1=step*f4(r,theta,p_r,p_theta);
            v1=step*f5(r,theta,p_r,p_theta);

            k2=step*f1(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2);
            m2=step*f2(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2);
            n2=step*f3(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2);
            s2=step*f4(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2);
            v2=step*f5(r+k1/2, theta+m1/2, p_r+n1/2, p_theta+s1/2);

            k3=step*f1(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2);
            m3=step*f2(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2);
            n3=step*f3(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2);
            s3=step*f4(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2);
            v3=step*f5(r+k2/2, theta+m2/2, p_r+n2/2, p_theta+s2/2);

            k4=step*f1(r+k3, theta+m3, p_r+n3, p_theta+s3);
            m4=step*f2(r+k3, theta+m3, p_r+n3, p_theta+s3);
            n4=step*f3(r+k3, theta+m3, p_r+n3, p_theta+s3);
            s4=step*f4(r+k3, theta+m3, p_r+n3, p_theta+s3);
            v4=step*f5(r+k3, theta+m3, p_r+n3, p_theta+s3);

            r=r+(k1+(2*k2)+(2*k3)+k4)/6;
            theta=theta+(m1+(2*m2)+(2*m3)+m4)/6;
            p_r=p_r+(n1+(2*n2)+(2*n3)+n4)/6;
            p_theta=p_theta+(s1+(2*s2)+(2*s3)+s4)/6;
            phi=phi+(v1+(2*v2)+(2*v3)+v4)/6;

            k = k+1;
            disp(r)

            x = [r theta p_r p_theta phi];
            curve(k,:)=x;



            %x=[r theta p_r p_theta];
            %disp(x);
            %disp(r);
            %disp(theta);
            %disp(["k=",k]);


        end

        %Transform to euclidean coordinates and plot3
        [n,m]=size(curve);
        A=a*ones(n,1);

        PHIZ=ones(n,1);

        Boyer2Cart=@(r,theta,phi)[sqrt(r.^2+A.^2).*sin(theta).*cos(phi),...
            sqrt(r.^2+A.^2).*sin(theta).*sin(phi),...
            r.*cos(theta)];

        cart=Boyer2Cart(curve(:,1),curve(:,2),curve(:,5));
%         disp([i,j]);
%         pause;

%         formatSpec = '%4.2f \t %4.2f \t %4.2f\n';
%         filename = sprintf('plot%d.txt', fig);
%         fileID = fopen(filename,'w');
%         fprintf(fileID,formatSpec,cart(:,1),cart(:,2),cart(:,3));  % The format string is applied to each element of a
%         fclose(fileID);
%         fig = fig+1


        plot3(cart(:,1),cart(:,2),cart(:,3))
        hold on


%         filename = sprintf('plot%d.png', fig);
%         saveas(figi, filename);
%         fig = fig+1;
         hold on;


   end
end

%close (hbar);
% if plot_type == 1
%     hold off
% end

%plot_r = curve(1,:)
%plot_theta = curve(2,:)
%polar(plot_theta,plot_r)
%A = [plot_theta',plot_r']

%disp(r);
%disp(theta);







