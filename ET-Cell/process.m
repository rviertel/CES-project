load ET.dat
load current.dat
figure(1)
plot(ET(:,1),ET(:,end)) %tau
title('tau')
figure(2)
plot(ET(:,1),ET(:,2))   %V
title('voltage')
figure(3)
plot(ET(:,1),ET(:,end-1))  %n
title('n')
