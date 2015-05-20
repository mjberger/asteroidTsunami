% plot gauge info in matlab.
% assume ngauges already set

ngauges = 7;
load _output/fort.gauge
figno =  64;  % starting figure number
clear legTrm;  % need new size to exactly correspond to number of gauges plotted, so reset
numPlots = 0;
legTrm = cellstr(num2str(zeros(ngauges,1)));

figure(figno); hold off
figure(figno+1); hold off
figure(figno+2); hold off
colors = ['b','g','r','k','c','m','k'];

for i= 1 : ngauges 
  numPlots = numPlots+1;
  index = find(fort(:,1) == i);
  gaugedata = fort(index,:);
  dp = gaugedata(:,8);
  eta = gaugedata(:,7);
  time = gaugedata(:,3);

  uvel = gaugedata(:,5)./gaugedata(:,4);
  vvel = gaugedata(:,6)./gaugedata(:,4);
  q = sqrt(uvel.*uvel+vvel.*vvel);

  numColors = 7;
  colIndex = mod(i,numColors);
  if (colIndex == 0) colIndex = numColors;end;
  col = colors(colIndex);
  
  figure(figno);plot(time,eta,col);
  hold on
  figure(figno+1);plot(time,dp,col);
  hold on
  figure(figno+2);plot(time,q,col);
  hold on  
  
  nextTerm = cellstr(sprintf(['gauge ',num2str(i)]));
  legTrm(numPlots) =  nextTerm;
end

figure(figno)
xlabel('time (sec)');
ylabel('eta (m)');
title('wave height at gauges');
%axis([0 .25 -3.0 .8])
grid on
legend(legTrm(1:numPlots))

figure(figno+1)
xlabel('time (sec)');
ylabel('dp/p');
title('pressure ratio at gauges');
axis([min(time) 1.2*max(time) -.11 .11])
grid on
legend(legTrm(1:numPlots))



figure(figno+2)
xlabel('time (sec)');
ylabel('velocity');
title('velocity at gauges');
%axis([0 max(time) -.002 .012])
grid on
legend(legTrm(1:numPlots))

