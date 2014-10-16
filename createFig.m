
% Phase 1 compare gm

% correlation .6
data = [ ...
    0.0180    0.0148    0.0062    0.0044
    0.0075    0.0048    0.0020    0.0033
    0.0075    0.0073    0.0045    0.0033
    0.0230    0.0118    0.0078    0.0070
    0.0390    0.0143    0.0078    0.0063
    0.2410    0.0125    0.0058    0.0047
    0.0420    0.0273    0.0168    0.0108
];
x = 2*[50   100   150   200];

h = figure
plot(x, data)
title('Comparison of GM methods')
xlabel('Number of Vertices')
ylabel('Accuracy')
legend('FAQ', 'Chance', 'Umeyama', 'Rank', 'QCV', 'PATH', 'GLAG');

saveas(gcf, 'Phase1_0.6.fig')
%savefig(h, 'Phase1_0.6.fig')



% correlation .9 %missing glag
data = [ ...
    0.1855    0.1180    0.0583    0.0143
    0.0075    0.0048    0.0020    0.0033
    0.0395    0.0160    0.0083    0.0074
    0.0515    0.0290    0.0167    0.0135
    0.9665    0.1760    0.0178    0.0138
    1.0000    0.3152    0.3583    0.0632
    0.8250    0.5500    0.3672    0.2719
];
x = 2*[50   100   150   200];

h = figure
plot(x, data)
title('Comparison of GM methods')
xlabel('Number of Vertices')
ylabel('Accuracy')
legend('FAQ', 'Chance', 'Umeyama', 'Rank', 'QCV', 'PATH', 'GLAG');

saveas(gcf, 'Phase1_0.9.fig')
%savefig(h, 'Phase1_0.9.fig')


% phase 2
% correlation .6
data = [ ...
    0.4312    0.4155    0.4603    0.5308    0.5872
    0.0097    0.0057    0.0043    0.0042    0.0026
    0.0090    0.0041    0.0034    0.0029    0.0013
    0.0099    0.0043    0.0039    0.0028    0.0021
    0.0135    0.0086    0.0050    0.0043    0.0036
    0.0163    0.0127    0.0081    0.0063    0.0038
];
x = [100, 200, 300, 400, 500];

h = figure
plot(x, data)
title('Comparison of LSGM methods')
xlabel('Max Cluster Size')
ylabel('Accuracy')
legend('SGM','FAQ', 'Chance', 'Umeyama', 'Rank', 'QCV');

saveas(gcf, 'Phase2_0.6.fig')
%savefig(h, 'Phase2_0.6.fig')


% correlation .9
data = [ ...
    0.6796    0.6129    0.6834    0.6774    0.7073
    0.0478    0.0173    0.0092    0.0075    0.0047
    0.0126    0.0046    0.0035    0.0025    0.0023
    0.0170    0.0092    0.0052    0.0048    0.0037
    0.0285    0.0153    0.0119    0.0091    0.0059
    0.1575    0.1141    0.0780    0.0592    0.0161
];
x = [100, 200, 300, 400, 500];

h = figure
plot(x, data)
title('Comparison of LSGM methods')
xlabel('Max Cluster Size')
ylabel('Accuracy')
legend('SGM','FAQ', 'Chance', 'Umeyama', 'Rank', 'QCV');

saveas(gcf, 'Phase2_0.9.fig')
%savefig(h, 'Phase2_0.9.fig')

% EXCLUDING LSGM in plot (too high)
%
% correlation .6
data = [ ...
%    0.4312    0.4155    0.4603    0.5308    0.5872
    0.0097    0.0057    0.0043    0.0042    0.0026
    0.0090    0.0041    0.0034    0.0029    0.0013
    0.0099    0.0043    0.0039    0.0028    0.0021
    0.0135    0.0086    0.0050    0.0043    0.0036
    0.0163    0.0127    0.0081    0.0063    0.0038
];
x = [100, 200, 300, 400, 500];

h = figure
plot(x, data)
title('Comparison of LSGM methods')
xlabel('Max Cluster Size')
ylabel('Accuracy')
legend('FAQ', 'Chance', 'Umeyama', 'Rank', 'QCV');

saveas(gcf, 'Phase2_0.6_no_lsgm.fig')
%savefig(h, 'Phase2_0.6.fig')


% correlation .9
data = [ ...
%    0.6796    0.6129    0.6834    0.6774    0.7073
    0.0478    0.0173    0.0092    0.0075    0.0047
    0.0126    0.0046    0.0035    0.0025    0.0023
    0.0170    0.0092    0.0052    0.0048    0.0037
    0.0285    0.0153    0.0119    0.0091    0.0059
    0.1575    0.1141    0.0780    0.0592    0.0161
];
x = [100, 200, 300, 400, 500];

h = figure
plot(x, data)
title('Comparison of LSGM methods')
xlabel('Max Cluster Size')
ylabel('Accuracy')
legend('FAQ', 'Chance', 'Umeyama', 'Rank', 'QCV');

saveas(gcf, 'Phase2_0.9_no_lsgm.fig')
%savefig(h, 'Phase2_0.9.fig')






