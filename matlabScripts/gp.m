maxAmp = 800;
width=90;
thick=12;
speed=.3915;
p_t = thick*2;
c=width/2.35482;
figure; 
%hold off
clear press
clear d

%time = 50;

for time = 1:1:200
%for time=65:65   
t = time*speed;
itno = 0;

  for i=0:.125:100
      itno = itno + 1;
      dist = i;
      d(itno) = i;
      g = maxAmp * exp(-15*(dist/c)^2)/100.;
      tail = 1./(1.+(.05*dist)^2);
      amplitude=g+tail;
      press(itno) = 0.;

      if (dist<t)
         press(itno) = amplitude* exp(-5.0*(t-dist)/p_t) *     ...
         (1.0 - 5.0*(t-dist)/p_t);

      end
   end
   plot(d,press);hold on;
     
end
grid on;

