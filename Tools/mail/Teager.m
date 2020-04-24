function en=Teager(sig)


Lsig  = length(sig);

en = sig(2:Lsig-1).^2 - sig(1:Lsig-2).*sig(3:Lsig);

end