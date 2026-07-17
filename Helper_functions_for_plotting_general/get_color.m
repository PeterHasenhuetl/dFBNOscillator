function color = get_color
% Specifies a set of colors as RGB codes for plotting.
% Part of a set of functions to simplify the making of figures in MATLAB
% (R2024a).
% Code written by Peter Hasenhuetl.


% First, individual colors in RGB
color.grey_blue = [0.5 0.6 0.7];
color.pastel_red = [1 0.6 0.7];
color.dark_blue = [0.095 0.3 0.5];
color.salmon = [1 0.35 0.37];
color.red = [0.4 0 0];

color.sleep_deprived = [240, 138, 162]./255;
color.rested = ([146, 198, 234])./255;

color.blue = [0 0 0.4];
color.daylight = [1 1 0.37];
color.gray = [0.45 0.45 0.45];
color.purple1 = [85 0 115]/255;
color.purple2 = [180  53 180]/255;
color.light_lavender = [202 158 174]/255;
color.dark_lavender = [157 111 127]/255;
color.light_yellow = [255 255 95]/255;
color.yellow = [252 197 78]/255;
color.dark_yellow = [255 165 0]/255;
color.orange = [255 100 0]/255;
color.magenta = [215 0 155]/255;
color.grass_green = [157 170 127]/255;
color.light_green = [0 191.25 0]/255;
color.medium_green = [0 127.5 0]/255;
color.dark_green = [0 63.75 0]/255;
color.navy = [66 109 127]/255;
color.red = [215 0 0]/255;
color.maroon = [128 0 0]/255;
color.peach = [255 176 124]/255;
color.smooth_pink1 = [255 110 135]/255;
color.smooth_pink2 = [255 120 175]/255;
color.smooth_pink3 = [220 120 150]/255;

color.teal = [0 128 128]/255;
color.light_gray = [191.25 191.25 191.25]/255;
color.medium_gray = [127.5 127.5 127.5]/255;
color.dark_gray = [63.75 63.75 63.75]/255;

color.blue_extra1 = [0 171.25 190]/255;
color.blue_extra2 = [0 191.25 190]/255;
color.light_blue1 = [110 150 255]/255;
color.light_blue2 = [110 120 240]/255;
color.light_blue3 = [140 153 255]/255;

color.greenish = [0 153 153]/255;
color.dark_greenish = [0 89.25 89.25]/255;

color.pruple = [85/255 0 115/255];
color.color1 = [202/255 158/255 174/255];
color.color2 = [202/255 158/255 63/255];
color.color3 = [157/255 111/255 127/255];
color.color4 = [157/255 170/255 127/255];
color.color5 = [222/255 177/255 170/255];
color.color6 = [128/255 109/255 127/255];
color.color7 = [66/255 109/255 127/255];

upward_blue_map(:,1) = linspace(1,0,10000)';
upward_blue_map(:,2) = linspace(1,0,10000)';
upward_blue_map(:,3) = ones(10000,1)';
color.upward_blue_map = upward_blue_map;

upward_grey_blue_map(:,1) = linspace(0.5,0,10000)';
upward_grey_blue_map(:,2) = linspace(0.6,0,10000)';
upward_grey_blue_map(:,3) = (ones(10000,1)').*0.7;
color.upward_grey_blue_map = upward_grey_blue_map;



%% Color maps for heat maps


downward_red_map(:,1) = ones(10000,1)';
downward_red_map(:,2) = linspace(0,1,10000)';
downward_red_map(:,3) = linspace(0,1,10000)';
color.downward_red_map = downward_red_map;


downward_blue_map(:,3) = ones(10000,1)';
downward_blue_map(:,2) = linspace(0,1,10000)';
downward_blue_map(:,1) = linspace(0,1,10000)';
color.downward_blue_map = downward_blue_map;

downward_gb_map(:,1) = linspace(0.5,1,10000)';
downward_gb_map(:,2) = linspace(0.6,1,10000)';
downward_gb_map(:,3) = linspace(0.99,1,10000)';
color.downward_gb_map = downward_gb_map;

downward_gb_map(:,1) = linspace(0.2,1,10000)';
downward_gb_map(:,2) = linspace(0.3,1,10000)';
downward_gb_map(:,3) = linspace(0.85,1,10000)';
color.downward_gb_map = downward_gb_map;


downward_salmon_map(:,1) = ones(10000,1)';
downward_salmon_map(:,2) = linspace(0,0.75,10000)';
downward_salmon_map(:,3) = linspace(0,0.7,10000)';
color.downward_salmon_map = downward_salmon_map;


down_bla_g_m(:,1) = zeros(10000,1)';
down_bla_g_m(:,2) = linspace(0,0.6,10000)';
down_bla_g_m(:,3) = linspace(0,0.6,10000)';
color.down_bla_g_m = down_bla_g_m;

up_bla_g_m(:,1) = zeros(10000,1)';
up_bla_g_m(:,2) = linspace(0.6,0,10000)';
up_bla_g_m(:,3) = linspace(0.6,0,10000)';
color.up_bla_g_m = up_bla_g_m;

upward_greenish_map(:,1) = linspace(0,1,10000)';
upward_greenish_map(:,2) = linspace(0.6,1,10000)';
upward_greenish_map(:,3) = linspace(0.6,1,10000)';
color.upward_greenish_map = upward_greenish_map;

down_greenish_map(:,1) = linspace(1,0,10000)';
down_greenish_map(:,2) = linspace(1,0.6,10000)';
down_greenish_map(:,3) = linspace(1,0.6,10000)';
color.down_greenish_map = down_greenish_map;

upward_green_map(:,1) = linspace(1,0,100)';
upward_green_map(:,2) = ones(100,1)';
upward_green_map(:,3) = linspace(1,0,100)';
color.upward_green_map = upward_green_map;

upward_black_map(:,1) = linspace(1,0,10000)';
upward_black_map(:,2) = linspace(1,0,10000)';
upward_black_map(:,3) = linspace(1,0,10000)';
color.upward_black_map = upward_black_map;

downward_black_map(:,1) = linspace(0,1,10000)';
downward_black_map(:,2) = linspace(0,1,10000)';
downward_black_map(:,3) = linspace(0,1,10000)';
color.downward_black_map = downward_black_map;

downward_blue_black_map(:,1) = zeros(10000,1);
downward_blue_black_map(:,2) = zeros(10000,1);
downward_blue_black_map(:,3) = linspace(1,0,10000)';
color.downward_blue_black_map = downward_blue_black_map;

upward_blue_black_map(:,1) = zeros(10000,1);
upward_blue_black_map(:,2) = zeros(10000,1);
upward_blue_black_map(:,3) = linspace(0,1,10000)';
color.upward_blue_black_map = upward_blue_black_map;

upward_yellow_black_map(:,2) = [zeros(2000,1); linspace(0,1,18000)'];
upward_yellow_black_map(:,1) = [linspace(0,1,2000)'; ones(18000,1)];
upward_yellow_black_map(:,3) = [linspace(0,1,18000)'; ones(2000,1)];
color.mixed4 = upward_yellow_black_map;


upward_yellow_black_map(:,1) = [linspace(0,0.75,2000)'; ones(18000,1).*0.75];
upward_yellow_black_map(:,2) = [linspace(0,0.75,2000)'; ones(18000,1).*0.75];
upward_yellow_black_map(:,3) = linspace(0,1,20000)';
color.mixed4 = upward_yellow_black_map;

upward_red_map(:,1) = ones(10000,1)';
upward_red_map(:,2) = linspace(1,0,10000)';
upward_red_map(:,3) = linspace(1,0,10000)';
color.upward_red_map = upward_red_map;

downward_red_black_map(:,1) = zeros(10000,1);
downward_red_black_map(:,3) = zeros(10000,1);
downward_red_black_map(:,1) = linspace(1,0,10000)';
color.mixed5 = downward_red_black_map;

salmon_map(:,1) = [ones(10000,1); linspace(1,0.3,10000)'];
salmon_map(:,2) = [linspace(1,0.35,10000)'; linspace(0.35,0.3,10000)'];
salmon_map(:,3) = [linspace(1,0.5,10000)'; linspace(0.5,0.3,10000)'];

upward_dark_blue_map(:,1) = [linspace(1,0,10000)'; zeros(100,1)];
upward_dark_blue_map(:,2) = [linspace(1,0,10000)'; zeros(100,1)];
upward_dark_blue_map(:,3) = [ones(10000,1); linspace(1,0,100)'];
color.upward_dark_blue_map = upward_dark_blue_map;

downward_dark_blue_map(:,1) = [linspace(0,1,1000)'; ones(10000,1)];
downward_dark_blue_map(:,2) = [linspace(0,1,1000)'; ones(10000,1)];
downward_dark_blue_map(:,3) = [zeros(1000,1); linspace(0,1,10000)'];
color.downward_dark_blue_map = downward_dark_blue_map;


color.mixed = [downward_gb_map; salmon_map];
color.mixed1 = [upward_blue_map; downward_blue_black_map];
color.mixed2 = [down_greenish_map; up_bla_g_m];
color.mixed3 = [down_bla_g_m(3001:end,:); upward_greenish_map];
color.mixed6 = [upward_blue_black_map(501:end,:);downward_blue_map; upward_red_map; downward_red_black_map(1:9500,:)];
color.mixed6_1 = [upward_blue_black_map(1:end,:);downward_blue_map; upward_red_map; downward_red_black_map(1:10000,:)];





%%

length_blue = 1000;

short_downward_blue_map(:,3) = ones(length_blue,1)';
short_downward_blue_map(:,2) = linspace(0,1,length_blue)';
short_downward_blue_map(:,1) = linspace(0,1,length_blue)';

short_upward_blue_black_map(:,1) = zeros(length_blue,1);
short_upward_blue_black_map(:,2) = zeros(length_blue,1);
short_upward_blue_black_map(:,3) = linspace(0,1,length_blue)';

color.skewed_blue_red = [short_upward_blue_black_map(1:end,:);short_downward_blue_map; upward_red_map];

%%

down_purple_map(:,1) = linspace(0.3333,1,20000)';
down_purple_map(:,2) = linspace(0,1,20000)';
down_purple_map(:,3) = linspace(0.4510,1,20000)';

color.down_purple_map = down_purple_map;
color.mixed_green_purple = [down_purple_map; down_greenish_map; up_bla_g_m];

down_purple_map_short(:,1) = linspace(0.3333,1,5000)';
down_purple_map_short(:,2) = linspace(0,1,5000)';
down_purple_map_short(:,3) = linspace(0.4510,1,5000)';
color.mixed_green_purple_skewed = [down_purple_map_short; down_greenish_map; up_bla_g_m];

%%

a(:,1) = linspace(0.95,0,5000)';
a(:,2) = linspace(0.95,0,5000)';
a(:,3) = ones(5000,1);

downward_blue_black_m(:,1) = zeros(5000,1);
downward_blue_black_m(:,2) = zeros(5000,1);
downward_blue_black_m(:,3) = linspace(1,0,5000)';
color.downward_blue_black_m = downward_blue_black_m;

upward_blue_black_m(:,3) = zeros(5000,1);
upward_blue_black_m(:,2) = zeros(5000,1);
upward_blue_black_m(:,1) = linspace(0,1,5000)';
color.upward_blue_black_m = upward_blue_black_m;

downward_blue_m(:,1) = ones(5000,1)';
downward_blue_m(:,2) = linspace(0,1,5000)';
downward_blue_m(:,3) = linspace(0,1,5000)';
color.downward_blue_m = downward_blue_m;

color.map_dFB_circadian = [a; downward_blue_black_m; upward_blue_black_m; downward_blue_m];

color.map_sleep = [[1,1,1]; [0.85, 0.85, 0.85]; [0,0,0]; [0,0,0]];


%%


dm(:,1) = zeros(20000,1);
dm(:,2) = zeros(20000,1);
dm(:,3) = linspace(0,0.75,20000)';

um(:,1) = linspace(0,1,225000)';
um(:,2) = linspace(0,1,225000)';
um(:,3) = linspace(0.75,0,225000)';


um(:,1) = linspace(0,1,225000)';
um(:,2) = linspace(0,0.85,225000)';
um(:,3) = zeros(225000,1);

dm(:,1) = [linspace(0.6,0,10000)';zeros(10000,1)];
dm(:,2) = [linspace(0.6,0,10000)';zeros(10000,1)];
dm(:,3) = linspace(1,0,20000)';

color.mixed7 = [dm; um];
color.color2 = color.color2.*1.25;

color.color5 = color.dark_blue;
color.mixed10 = [linspace(color.color2(1),color.color7(1),10000)',...
    linspace(color.color2(2),color.color7(2),10000)', linspace(color.color2(3),color.color7(3),10000)'];

color.mixed10 = [linspace(color.dark_blue(1),color.salmon(1),10000)',...
    linspace(color.dark_blue(2),color.salmon(2),10000)', linspace(color.dark_blue(3),color.salmon(3),10000)'];


color.mixed11_1 = [linspace(color.color2(1),1,10000)',...
    linspace(color.color2(2),1,10000)', linspace(color.color2(3),1,10000)'];

color.mixed11_2 = [linspace(color.salmon(1),color.color2(1),10000)',...
    linspace(color.salmon(2),color.color2(2),10000)', linspace(color.salmon(3),color.color2(3),10000)'];

color.mixed11_3 = [linspace(color.dark_blue(1),color.salmon(1),10000)',...
    linspace(color.dark_blue(2),color.salmon(2),10000)', linspace(color.dark_blue(3),color.salmon(3),10000)'];

color.mixed11_4 = [linspace(0,color.dark_blue(1),10000)',...
    linspace(0,color.dark_blue(2),10000)', linspace(0,color.dark_blue(3),10000)'];

color.mixed10 = [color.mixed11_4; color.mixed11_3; color.mixed11_2; color.mixed11_1];


%%

d_black = linspace(1,0,10000)';
color.plusdF = color.mixed10;
color.minusdF = [d_black, d_black, d_black];


end
