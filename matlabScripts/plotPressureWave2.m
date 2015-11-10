%% plot pressure field model (widht and ampitude) 
%% for chelyabinsk and tunguska

%% dividing everything by 2 reduced the underpressure, but will still approach vacuum
%% when added to atmospheric pressure


tunguska = true;

dist_in_km = 0:.01:100;
if (tunguska)
   figno = 56;
   maxAmp = 800;  % 179; in percent
   width = 90;  %50;
   str = 'Tunguska; time (sec) ';
else
    figno = 45;
    maxAmp = 06.5;
    width=150;
    str = 'Chelyabinsk; time (sec) ';
end
thick = 5;
speed = .3915;  % in km/sec

%% pick  times to evaluate pressure wave
time = [50 100 150 200 250];  % in 100 sec it moves 39km
%time = [10 25 50  100 150 ];  % in 100 sec it moves 39km

figure(figno); hold off;

%for it=1:length(time)
for it=1:5:300
   %t = time(it) * speed;
   t = it*speed;
   
   p_t  =  thick*2.d0;
   c    =  width/2.35482d0;
   g    =  maxAmp .* exp(-dist_in_km.^2./(c*c));  % vector op for all dists
   computedOverPressure = g.*exp(-0.8d0*(t-dist_in_km)/p_t).*   ...
        (1.d0 - 1.1d0*(t-dist_in_km)/p_t);

  computedOverPressure = computedOverPressure/2.;  % for new sim with half max mpl
  % just take half of everything, so min is halved too. min went too negative

   index = find(dist_in_km >= t);
   computedOverPressure(index) = 0.d0;  % if wave hasn't reached there yet
   
   hold on
   %figure(figno+it);plot(dist_in_km,computedOverPressure);
   %timeString = sprintf([str,num2str(time(it))]);
   %legend(timeString);

   plot(dist_in_km,computedOverPressure);
   hold on;plot(-dist_in_km,computedOverPressure);

   xlabel('Dist (km)');
   ylabel('% Overpressure (atm)');
   grid on;
end

%legend('50 sec', '100 sec','150 sec', '200 sec', '250 sec');
%legend('10 sec','25 sec', '50 sec','75 sec','100 sec', '150 sec');
%axis([-100. 100. -4.1 10.])
