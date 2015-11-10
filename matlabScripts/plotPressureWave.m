%% plot pressure field model (widht and ampitude) 
%% for chelyabinsk and tunguska


%% way too large underpressure - will make atmospheric pressure go negative 


tunguska = true;

dist_in_km = 0:.1:100;
if (tunguska)
   figno = 55;
   maxAmp = 800;  % 179; in percent
   width = 90;  %50;
   str = 'Tunguska; time (sec) ';
else
    figno = 45;
    maxAmp = 06;
    width=150;
    str = 'Chelyabinsk; time (sec) ';
end
thick = 5;
speed = .3915;  % in km/sec

%% pick  times to evaluate pressure wave
time = [25 50 100 150 200];  % in 100 sec it moves 39km
%time = [10 25 50  100 150 ];  % in 100 sec it moves 39km

figure(figno); hold off;

%for it=1:length(time)
for it=1:5:200
    
   %t = time(it) * speed;
    t = it*speed;
   
   p_t  =  thick*2.d0;
   c    =  width/2.35482d0;
   g    =  maxAmp .* exp(-dist_in_km.^2./(c*c));  % vector op for all dists
   computedOverPressure = 2.d0*g.*exp(-0.8d0*(t-dist_in_km)/p_t).*   ...
        (0.5d0 - 1.1d0*(t-dist_in_km)/p_t);

   index = find(dist_in_km >= t);
   computedOverPressure(index) = 0.d0;  % if wave hasn't reached there yet
   
   hold on
   %figure(figno+it);plot(dist_in_km,computedOverPressure);
   %timeString = sprintf([str,num2str(time(it))]);
   %legend(timeString);

   plot(dist_in_km,computedOverPressure);

   xlabel('Dist (km)');
   ylabel('% Overpressure (atm)');
   grid on;
end

%legend('25 sec', '50 sec','100 sec', '150 sec', '200 sec');
%legend('10 sec','25 sec', '50 sec','75 sec','100 sec', '150 sec');

