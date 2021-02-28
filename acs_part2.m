
time=[0 15];
yo=[0 0 0 0];
%\\\\\\\\\\\\\\\\\\\ERWTHMA A \\\\\\\\\\\\\\\\\\\
[t,x]=ode45(@func1,time,yo);
plot(t,x(:,1),'c')
title('Position')
legend('x1')
figure
plot(t,x(:,2),'k')
title('Velocity')
legend('x2')
figure
plot(t, x(:,1) - 0.1 - 0.5*sin(t) , 'g' )
title('Error of position')
legend('m')
figure
plot(t,x(:,2) - 0.5*cos(t) ,'b')
title('Error of velocity')
legend('edot')

V=-mi/J.*x(:,4)-k/J.*(x(:,3)-x(:,1))+1/J*0.1 + 0.5*sin(t)- ko.*(x(:,1)-0.1-0.5*sin(t)) -k1.*(x(:,2)-0.5*cos(t))-k2.*(-(m*g*r)/I*cos(x(:,1))+k/I.*(x(:,3)-x(:,1))+0.5*sin(t))-k3.*((m*g*r)/I.*sin(x(:,1)).*x(:,2)+k/I.*(x(:,4)-x(:,2))+0.5*cos(t)) ;
 
figure
plot(t, V ,'r')
title('V')
legend('V')
 
de = interp1(t,x(:,1),10);
dy = interp1(t,0.1 + 0.5*sin(t),10);
error = abs(de - dy);
display(error);

%\\\\\\\\\\\\\\\\\\\ERWTHMA B \\\\\\\\\\\\\\ EKTIMHSEIS

[t,x]=ode45(@func2,time,yo);
plot(t,x(:,1),'c')
title('Position')
legend('x1')
figure
plot(t,x(:,2),'k')
title('Velocity')
legend('x2')
figure
plot(t, x(:,1) - 0.1 - 0.5*sin(t) , 'g' )
title('Error of position')
legend('m')
figure
plot(t,x(:,2) - 0.5*cos(t) ,'b')
title('Error of velocity')
legend('edot')

de = interp1(t,x(:,1),10);
dy = interp1(t,0.1 + 0.5*sin(t),10);
error = abs(de - dy);
display(error);

  I_hat=9;
  m_hat=23;
  r_hat=0.4;
  k_hat=2400;
  mi_hat=2.5;
  J_hat=9;
    
    
 V=-mi_hat/J_hat.*x(:,4)-k_hat/J_hat.*(x(:,3)-x(:,1))+1/J_hat*0.1 + 0.5*sin(t)- ko.*(x(:,1)-0.1-0.5*sin(t)) -k1.*(x(:,2)-0.5*cos(t))-k2.*(-(m_hat*g*r_hat)/I_hat*cos(x(:,1))+k_hat/I_hat.*(x(:,3)-x(:,1))+0.5*sin(t))-k3.*((m_hat*g*r_hat)/I_hat.*sin(x(:,1)).*x(:,2)+k/I_hat.*(x(:,4)-x(:,2))+0.5*cos(t)) ;
 
figure
plot(t, V ,'red')
title('V')
legend('V')

%\\\\\\\\\\\\\\ERWTHMA C      \\\\\\\\\\\\
options=odeset('Reltol',1e-6,'AbsTol',1e-3);
[t,x]=ode23s(@func3,time,yo,options);

de = interp1(t,x(:,1),10);
dy = interp1(t,0.1 + 0.5*sin(t),10);
error = abs(de - dy);
display(error);
    
  L = 2000 ; 
    
    e_dot3 = ( m*g*r.*sin(x(:,1)).*x(:,2) ) / I1 + (k1.*(x(:,4)-x(:,2)))/I1 + 0.5*cos(t);

    e_dot2 = ( -m*g*r.*cos(x(:,1)) +k1.*(x(:,3)-x(:,1)) )/I1 + 0.5*sin(t) ;

    e_dot1 = x(:,2) - 0.5*cos(t) ;

    e_dot0 = x(:,1) - 0.1 -0.5*sin(t) ;

    ueq = -((J1_hat*m_hat*g*r_hat)/k1_hat).*(cos(x(:,1)).*x(:,2).^2) -((J1_hat*m_hat*g*r_hat)/k1_hat).*(sin(x(:,1))).*((-m_hat*g*r_hat.*cos(x(:,1)) + (k1_hat/I1_hat).*(x(:,3)-x(:,1)))/I1_hat) +M_hat.*x(:,4) -k1_hat.*(x(:,3)-x(:,1)) - ( ( J1_hat / k1_hat ) * (m_hat * g * r_hat .* cos(x(:,1)) - k1_hat .* ( x(:,3) - x(:,1)) ) ) + ( (I1_hat*J1_hat)/k1_hat).*(-0.5*cos(t))  - ((J1_hat*m_hat*g*r_hat)/k1_hat).*( 3.*e_dot3.*L + 3.*e_dot2.*(L^2) +(L^3).*e_dot1) ;                         
     
    s = e_dot3+3.*L.*e_dot2 + 3.*(L^2).*e_dot1 + (L^3).*e_dot0 ; 
    
    A = ((J1*m*g*r)/k1) - ((J1_hat*m_hat*g_hat*r_hat)/k1_hat);
    B = ((J1*(m*g*r)^2)/(k1*I1)) - ((J1_hat*(m_hat*g_hat*r_hat)^2)/(k1_hat*I1_hat));
    C =  ((J1*m*g*r)/(I1^2)) - ((J1_hat*m_hat*g_hat*r_hat)/(I1_hat^2));
    D = J1 - J1_hat;
    E = k1 - k1_hat ; 
    F = ((I1*J1)/k1) - ((I1_hat*J1_hat)/k1_hat);
    G = M-M_hat ;
    
    c = 10;
    d = 0.00001;
    
    value = min(max(s/d, -1), 1);
   
    R =  abs(cos(x(:,1))).*abs(x(:,2).^2).*max(abs(A1)) - abs(sin(x(:,1))).*abs(cos(x(:,1))).*max(abs(A2)) + abs(sin(x(:,1))).*abs(x(:,3)-x(:,1))*max(abs(A3)) + abs(cos(x(:,1))).*max(abs(A1)) - abs(x(:,3)-x(:,1)).*max(abs(A4)) - abs(x(:,4)).*max(abs(A7)) + abs(x(:,3)-x(:,1)).*max(abs(A5)) -(abs(-0.5*cos(t)))*max(abs(A6)) + abs( 3.*e_dot3.*L + 3.*e_dot2.*(L^2) +(L^3).*e_dot1 ).* max(abs(A6)) + c    ;
    
    u = ueq - R.*value;

    figure
    plot(t,u,'r')
    title('Response u')
    legend('u')
    figure

    plot(t,s,'g')
    legend('s')
    title('Function s')
   
% SYNARTHSH TOU ERWTHMATOS A , SXETIKA ME TO ERROR 
    function xdot=func1(t,x)
    
    I=5.6;
    m=20;
    g=9.8;
    r=0.25;
    k=1900;
    mi=1;
    J=6.18;
    
    z(1)=x(1);
    z(2)=x(2);
    z(3)=-(m*g*r)/I*cos(x(1))+k/I*(x(3)-x(1));
    z(4)=(m*g*r)/I*sin(x(1))*x(2)+k/I*(x(4)-x(2));
    %kerdi
    ko=10000;
    k1=4000;
    k2=600;
    k3=40;
    
    xdot(1)=x(2);
    xdot(2)=-(m*g*r)/I*cos(x(1))+k/I*(x(3)-x(1));
    xdot(3)=x(4);
    xdot(4)=-mi/J*x(4)-k/J*(x(3)-x(1))+1/J*0.1 + 0.5*sin(t)- ko*(z(1)-0.1-0.5*sin(t)) -k1*(z(2)-0.5*cos(t))-k2*(z(3)+0.5*sin(t))-k3*(z(4)+0.5*cos(t)) ;
    
    xdot=xdot';
    end
    % SYNARTHSH TOU ERWTHMATOS B POY PERIEXEI KAI TIS EKTIMHSEIS
 function xdot=func2(t,x)
    I=5.6;
    m=20;
    g=9.8;
    r=0.25;
    k=1900;
    mi=1;
    J=6.18;
    
    %ektimiseis
    I_hat=9;
    m_hat=23;
    r_hat=0.4;
    k_hat=2400;
    mi_hat=2.5;
    J_hat=9;
    
    z(1)=x(1);
    z(2)=x(2);
    z(3)=-(m_hat*g*r_hat)/I_hat*cos(x(1))+k/I*(x(3)-x(1));
    z(4)=(m_hat*g*r_hat)/I_hat*sin(x(1))*x(2)+k/I*(x(4)-x(2));
    %kerdi
    ko=10000;
    k1=4000;
    k2=600;
    k3=40;
    
    xdot(1)=x(2);
    xdot(2)=-(m*g*r)/I*cos(x(1))+k/I*(x(3)-x(1));
    xdot(3)=x(4);
    xdot(4)=-mi_hat/J_hat*x(4)-k_hat/J_hat*(x(3)-x(1))+1/J_hat*0.1 + 0.5*sin(t)- ko*(z(1)-0.1-0.5*sin(t)) -k1*(z(2)-0.5*cos(t))-k2*(z(3)+0.5*sin(t))-k3*(z(4)+0.5*cos(t)) ;
    
    xdot=xdot';
 end
% SYNARTHSH TOY ERWTHMATOS C POY PERIEXEI KAI TON ELEGXO OLISTHISIS
function xdot=func3(t,x)
    
    I1=5.6;
    m=20;
    g=9.8;
    r=0.25;
    k1=1900;
    M=1;
    J1=6.18;
 %ektimiseis
    I1_hat=4;
    m_hat=21;
    g_hat = 9.8;
    r_hat=0.32 ;
    k1_hat= 2000 ;
    M_hat = 0.9;
    J1_hat = 7;
    
    L = 500 ; 
    %paragwgoi sfalmatos
    e_dot3 = ( m*g*r*sin(x(1))*x(2) ) / I1 + (k1*(x(4)-x(2)))/I1 + 0.5*cos(t);
    e_dot2 = ( -m*g*r*cos(x(1)) +k1*(x(3)-x(1)) )/I1 + 0.5*sin(t) ;
    e_dot1 = x(2) - 0.5*cos(t) ;
    e_dot0 = x(1) - 0.1 -0.5*sin(t) ;

    ueq = -((J1_hat*m_hat*g*r_hat)/k1_hat)*(cos(x(1))*x(2)^2) -((J1_hat*m_hat*g*r_hat)/k1_hat)*(sin(x(1)))*((-m_hat*g*r_hat*cos(x(1)) + (k1_hat/I1_hat)*(x(3)-x(1)))/I1_hat) +M_hat*x(4) -k1_hat*(x(3)-x(1)) - ( ( J1_hat / k1_hat ) * (m_hat * g * r_hat * cos(x(1)) - k1_hat * ( x(3) - x(1)) ) ) + ( (I1_hat*J1_hat)/k1_hat)*(-0.5*cos(t))  - ((J1_hat*m_hat*g*r_hat)/k1_hat)*( 3*e_dot3*L + 3*e_dot2*(L^2) +(L^3)*e_dot1) ;                         
     
    s = e_dot3+3*L*e_dot2 + 3*(L^2)*e_dot1 + (L^3)*e_dot0 ; 
    %metavlites gia to R
    A = ((J1*m*g*r)/k1) - ((J1_hat*m_hat*g_hat*r_hat)/k1_hat);
    B = ((J1*(m*g*r)^2)/(k1*I1)) - ((J1_hat*(m_hat*g_hat*r_hat)^2)/(k1_hat*I1_hat));
    C =  ((J1*m*g*r)/(I1^2)) - ((J1_hat*m_hat*g_hat*r_hat)/(I1_hat^2));
    D = J1 - J1_hat;
    E = k1 - k1_hat ; 
    F = ((I1*J1)/k1) - ((I1_hat*J1_hat)/k1_hat);
    G = M-M_hat ;
    
    c = 10;
    d = 0.001;
    value = min(max(s/d, -1), 1);
   
    R =  abs(cos(x(1)))*abs(x(2)^2)*max(abs(A)) - abs(sin(x(1)))*abs(cos(x(1)))*max(abs(B)) + abs(sin(x(1)))*abs(x(3)-x(1))*max(abs(C)) + abs(cos(x(1)))*max(abs(A)) - abs(x(3)-x(1))*max(abs(D)) - abs(x(4))*max(abs(G)) + abs(x(3)-x(1))*max(abs(E)) -(abs(-0.5*cos(t)))*max(abs(F)) + abs( 3*e_dot3*L + 3*e_dot2*(L^2) +(L^3)*e_dot1 )* max(abs(F)) + c    ;
    u = ueq - R*value;
    
    xdot(1)=x(2);
    xdot(2)=(-m*g*r*cos(x(1)) + k1*(x(3)-x(1)))/I1;
    xdot(3)=x(4);
    xdot(4)= (-k1*(x(3)-x(1)) -M*x(4))/J1 + u/J1;
    
    xdot=xdot';
end
