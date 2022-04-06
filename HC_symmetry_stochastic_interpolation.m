function sdot= HC_symmetry_stochastic_interpolation(s,t,m1f,m1r,m2,kt,kf,kr,cr,cf,J,g,hf,hr,ttf,ttr,lr,lf,L)
Front_exct = interp1(ttf,hf,t,'spline');
Rear_exct = interp1(ttf,hr,t,'spline');
sdot = [s(2);
    (-kf*(s(1)-s(5)-lf*s(7))-kt*(s(1)-Front_exct)-cf*(s(2)-s(6)-lf*s(8))-m1f*g)/m1f; 
      s(4);
    (-kr*(s(3)-s(5)+lr*s(7))-kt*(s(3)-Rear_exct)-cr*(s(4)+lr*s(8)-s(6))-m1r*g)/m1r; 
      s(6);
    (-kr*(s(5)-lr*s(7)-s(3))-kf*(s(5)-s(1)+lf*s(7))-cr*(s(6)-lr*s(8)-s(4))-cf*(s(6)-s(2)+lf*s(8))-m2*g)/m2;
      s(8);
    (+lr*kr*(s(5)-lr*s(7)-s(3))+lr*cr*(s(6)-lr*s(8)-s(4))-lf*kf*(s(5)+lf*s(7)-s(1))-lf*cf*(s(6)+lf*s(8)-s(2)))/J];

end
