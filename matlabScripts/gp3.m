% trying to fiddle with new mega case to erduce underpressure;

maxAmp1 = 8000;
maxAmp2 =  200;
maxAmp3 =  200;
width=40;  % width of envelope

thick=5;  % of waveform itself
speed=.3915;
p_t = thick*2;
c=width/2.35482;
figure; 
%hold off
clear press
clear d

for time = 1:5:100
      t = time*speed;
      dist = 0:.05:100;
      g = 0 + maxAmp1 * exp(-4.5*(dist./c).^2);
      g = g + maxAmp2 * exp(-2.0*(dist./(5.5d0*c)).^2);
      g = g + maxAmp3 * exp(-2.0*(dist./(3.d0*c)).^4);
      press = 0. * dist;

      index = find(dist<t);

      press(index) = max(.5*g(index).*exp(-3.8d0*(t-dist(index))/p_t).*(2.0d0 - ...
      1.2d0*(t-dist(index))/p_t),-50);

   plot(dist,press);hold on;
     
end
grid on;
xlabel('km')
ylabel('overpressure')

