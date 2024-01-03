chamber_diameter = [0.176,0.066,0.1,0.273];
chamber_length =  [0.342,.2,.356,.323];
throat_diamter = [.101,0.023,0.058,0.198];

chamber_area = pi.*(chamber_diameter./2).^2;
chamber_volume = chamber_area.*chamber_length;

L_star = chamber_volume ./ (pi*(throat_diamter./2).^2);

avg_L_Star = mean(L_star);

our_chamber_v = 0.00500*avg_L_Star;

chamber_area = 1.1*0.00500;


chamber_length = our_chamber_v /chamber_area\

%%

flowisentropic(1.23986,)