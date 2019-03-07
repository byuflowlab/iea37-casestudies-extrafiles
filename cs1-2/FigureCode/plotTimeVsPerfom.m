% Nick Baker
% IEA37 Time vs Performance Plot
% 4 Dec 2018

% Directory for our figure
figuresdir = '/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/Figures/';


% Non-normalized data
% Data
% Participant order: (1)PJ, (2)Tim, (3)Erik full, (4)Jared, (5)Nick,
% (6)Sanchez, (7)Wiley, (8)Ning, (9)Erik partial, (10)Prakash
AEP = [1476689.663, 1506388.415, 1455075.608, 1513311.194, 1336164.550, 1364943.008, 1332883.433, 1445967.377, 1422268.714, 1480850.976];
AEP = AEP ./ 1e6;
Walltime = [288, 43.2, 10000, 5299.01, 57.6, 318, 1262.4, 90, 1000, 193];

% Plt
figure(1)
hold on
for i = 1:10
    p(i) = scatter(Walltime(i), AEP(i),'^', 'filled');
end
hold off
%legend(p,'PJ', 'Tim','Erik (full)','Thomas','Nick','Sanchez','Willey','Ning','Erick (simp)', 'Prakash', 'location', 'NorthWest')
legend(p,'SNOPT', 'Preconditioned Sequential Quadratic Programming','Full Pseudo-Gradient Approach','SNOPT+WEC','FMinCon()','Simple Particle Swarm Optimisation','Basic Genetic Algorithm','SNOPT','Simple Pseudo-Gradient Approach', 'Multistart Interior-Point', 'location', 'NorthWest')
xlabel('Reprted Wall Time (s)')
ylabel('AEP (MWh)')

% Normalized Graph
NormAEPvTime = AEP ./ Walltime;

figure(2)
hold on
for ParNum = 1:10
    q(i) = scatter(ParNum, NormAEPvTime(ParNum), '^', 'filled');
end
legend(q,'SNOPT', 'Preconditioned Sequential Quadratic Programming','Full Pseudo-Gradient Approach','SNOPT+WEC','FMinCon()','Simple Particle Swarm Optimisation','Basic Genetic Algorithm','SNOPT','Simple Pseudo-Gradient Approach', 'Multistart Interior-Point', 'location', 'NorthEast')
hold off
ylabel('MWH/s')
xlabel('sub#')





% Save the graph as a pdf
%saveas(gcf,strcat(figuresdir, fig_name))