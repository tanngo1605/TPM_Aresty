function dragCoefXY = getDragCoefXY(r, z, sita0)
ratio = r/z;
denom = 1 - (9/16)*ratio + (1/8)*ratio^3 - (45/256)*ratio^4 - (1/16)*ratio^5;
dragCoefXY = sita0 / denom;