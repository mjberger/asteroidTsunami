% extract line at given time from hovSlice.dat (also named all.dat in
% tecplot macro). plot 
%wave height  and overpressure in line plot

% must specify this to use script
%timeWanted =    53.789020000000001
timeWanted =    37.126910000000002
  
  
load ~/slides/tsunami/all.dat
x = all(:,1);
time=all(:,2);
eta = all(:,3);
overpress = all(:,4);  
%since didnt yet impleement outting of overpressure ratio
% convert now
overpress = overpress/101300 - 1;


% extract time
index = find(time==timeWanted);

hold off
plot(x(index),eta(index),'-+');
hold on
plot(x(index),overpress(index),'-x');

legend('height (m)','overpressure (atm)');
grid on;
