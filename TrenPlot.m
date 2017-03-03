% Plot Test



Tsolves = [24.1905089370000,26.7739301100000,29.1709141970000,31.2825510730000,37.7033725840000,31.0272098740000,33.2449475140000,32.8753750020000,32.9848359520000,35.0443795110000,37.2758347590000,44.2212441950000,72.1227674810000];


NumITS = [1954,2184,2355,2458,2508,2532,2554,2616,2753,2924,3107,3517,5805] ;


ReNums = [1 , 100, 200,300,400,500,600,700,800,900,1000,1100,1200] ;



figure(1)



plot(ReNums,NumITS,'r','LineWidth',2)

title('Itterations to reach steady state')
xlabel('Reynolds number')
ylabel('Number of Itterations')




figure(2)

plot(ReNums,Tsolves,'b','LineWidth',2)
title('Time taken to reach steady state')
xlabel('Reynolds number')
ylabel('Time')