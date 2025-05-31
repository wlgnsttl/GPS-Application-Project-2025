function [vec_sat_p, vec_sat_v] = getSatVel_diff(eph, ieph, tc, STT)
[vec_sat_p, ~] = getSatPos_lab(eph, ieph, tc);
[vec_sat_p_a, ~] = getSatPos_lab(eph, ieph, tc+1e-5);
vec_sat_p = RotSatPos(vec_sat_p, STT);
vec_sat_p_a = RotSatPos(vec_sat_p_a, STT);
vec_sat_p = vec_sat_p';
vec_sat_p_a = vec_sat_p_a';

vec_sat_v = (vec_sat_p_a-vec_sat_p)/1e-5;
end