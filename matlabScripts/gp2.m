% trying to fiddle with new mega case to erduce underpressure;

maxAmp = 8400;
width=45;  % width of envelope

thick=6;  % of waveform itself
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
      g = maxAmp * exp(-2*(dist./c).^2);
      press = 0. * dist;

      index = find(dist<t);

      press(index) = max(.5*g(index).*exp(-3.8d0*(t-dist(index))/p_t).*(2.0d0 - ...
      1.2d0*(t-dist(index))/p_t),-50);

   plot(dist,press);hold on;
     
end
grid on;
xlabel('km')
ylabel('overpressure')

