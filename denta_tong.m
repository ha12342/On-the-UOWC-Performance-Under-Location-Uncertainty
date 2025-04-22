function denta = denta_tong(x,y,z,uu_3,m)

denta = 2*(-uu_3)^(-(x+2/(m+1)))*(igamma(x+2/(m+1),-2*y/(m+1)) - igamma(x+2/(m+1),-2*z/(m+1)));

end