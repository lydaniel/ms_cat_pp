% ms
ms = [...
1, 0.0349225118603, 0.0073878456924; ...
2, 0.082580124979, 0.0034044075413; ...
3, 0.136459018478, 0.0014940968606; ...
4, 0.1918149799, 0.000790946587731; ...
5, 0.249267562522, 0.00090003437191; ...
6, 0.301242621425, 0.0009595452592; ...
7, 0.351176514524, 0.000984871137004; ...
8, 0.397737946578, 0.00098438740899; ...
9, 0.441423121855, 0.00094173858166; ...
10, 0.481720458105, 0.000883891183818; ...
11, 0.515308715304, 0.000820554240212; ...
12, 0.546415801165, 0.000731031886966; ...
13, 0.572758041753, 0.000653042281458; ...
14, 0.595351346135, 0.000560123169278; ...
15, 0.614256320747, 0.000485174339222; ...
16, 0.630769285329, 0.000393421258742; ...
17, 0.643685939528, 0.00032588585673; ...
18, 0.65441236102, 0.000263466985891; ...
19, 0.663244544492, 0.000206171700347; ...
20, 0.670623910724, 0.000160879021244; ...
21, 0.675942646049, 0.000122772736882; ...
22, 0.680339127958, 8.83775011828e-05; ...
23, 0.683443781653, 6.74642882698e-05; ...
24, 0.685949150954, 4.86094872477e-05; ...
25, 0.687906494867, 3.46418481779e-05; ...
26, 0.689316351861, 2.44792584073e-05; ...
27, 0.690346073227, 1.65554632988e-05; ...
28, 0.691157323714, 1.11026644765e-05; ...
29, 0.691715792493, 7.32335391029e-06; ...
30, 0.692129072059, 4.62614993687e-06; ...
31, 0.692431982947, 6.59350781077e-06; ...
32, 0.692647020515, 1.22271772874e-05; ...
33, 0.69279981269, 2.44929191806e-05; ...
34, 0.692907592865, 0; ...
];

%total
total = [ ...
1, 0.166822279603, 0.0280194533081; ...
2, 0.318905027782, 0.00724495489652; ...
3, 0.438640939844, 0.00186351522464; ...
4, 0.526248422511, 0.00107348978962; ...
5, 0.586360091044, 0.000880343980738; ...
6, 0.627183601949, 0.000680618551317; ...
7, 0.65280569657, 0.000526301144969; ...
8, 0.669628489427, 0.00037525591204; ...
9, 0.679691579228, 0.000266095203001; ...
10, 0.68571997977, 0.000185303506417; ...
11, 0.689158511353, 0.00012017859745; ...
12, 0.691030463497, 7.77836258925e-05; ...
13, 0.692027173178, 5.2026387402e-05; ...
14, 0.692572155432, 3.48033221459e-05; ...
15, 0.692854504264, 2.21559339341e-05; ...
16, 0.69300030037, 1.06532517389e-05; ...
17, 0.69309628409, 3.27854169655e-06; ...
18, 0.693123432322, 1.7692342858e-06; ...
19, 0.693135258123, 1.04038569633e-06; ...
20, 0.693141956062, 5.93421576093e-07; ...
21, 0.693145076195, 3.44728412576e-07; ...
22, 0.69314603968, 1.63121256544e-07; ...
23, 0.693146834071, 4.6552966297e-08; ...
24, 0.693147032719, 3.13939229497e-08; ...
25, 0.693147136312, 4.92168010051e-09; ...
26, 0.693147164066, 3.31169808476e-09; ...
27, 0.69314717499, 9.76618801302e-10; ...
28, 0.693147178177, 4.81486617669e-10; ...
29, 0.69314717983, 9.09156984193e-11; ...
30, 0.693147180336, 2.85092298826e-11; ...
31, 0.693147180493, 1.79339525032e-11; ...
32, 0.693147180539, 2.73754496744e-12; ...
33, 0.693147180553, 2.11648318143e-12; ...
34, 0.693147180558, 3.81095782375e-13; ...
35, 0.693147180559, 5.90728268602e-14; ...
36, 0.693147180559, 7.74284633617e-15; ...
37, 0.69314718056, 3.06920512868e-15; ...
38, 0.69314718056, 1.97319157609e-15; ...
39, 0.69314718056, 1.94944579505e-15; ...
40, 0.69314718056, 1.95234247267e-15; ...
41, 0.69314718056, 1.95367848332e-15; ...
42, 0.69314718056, 1.07379551428e-15; ...
43, 0.69314718056, 4.02454636757e-16; ...
44, 0.69314718056, 1.17160692928e-16; ...
45, 0.69314718056, 0; ...
];

hold on
errorbar(total(:,1), total(:,2), total(:,3), 'r');
errorbar(ms(:,1), ms(:,2), ms(:,3), 'b');

title('Empirical information gain');
ylabel('Information gain [nats] \rightarrow');
xlabel('Number of participants \rightarrow');
set(gcf, 'Position', [100, 100, 400, 300]);
%legend(gcf, {'Optimal experiment'; 'Medin and Schaffer 5-4'});
legend('Optimal experiment', 'Medin and Schaffer 5-4', 'Location', 'southeast');
axis([0 40 0 0.7]);
box on;

