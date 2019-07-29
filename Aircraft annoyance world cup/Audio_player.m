reference_nb = 9

load awc.mat
%aircraft data set selection based on reference number
signallst = {signal1,signal2,signal3,signal4,signal5,signal6,signal7,signal8,signal9,signal10,signal11,signal12,signal13};
signal = signallst{1,reference_nb};

sound(signal,40000)

