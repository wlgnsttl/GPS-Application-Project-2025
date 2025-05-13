function [fixed_xyz] = RotSatPos(xyz, STT)
w = 7.2921151467e-5;
tau=STT;
Euler_matrix = [cos(w*tau), sin(w*tau), 0;
                -sin(w*tau), cos(w*tau), 0;
                0, 0, 1];
rotation_fixed = Euler_matrix*xyz';
fixed_xyz = rotation_fixed';