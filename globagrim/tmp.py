# Case 2: Stream over an isolated montain (Williamson test similation)
def montain_flow():
    u0=20.
    global_const.phis0=9.81*2000.0
    flamc=3.*global_const.pi/2.
    phic=global_const.pi/6.
    distm=global_const.pi/9.
    for k in range (0,global_const.NK)
        for j in range (0,global_const.NJ)
            dist=sqrt((flam(j)-flamc)**2+(phi(k)-phic)**2)
            if(dist.lt.distm):
                global_const.phis(j,k)=global_const.phis0*(1.-dist/distm)
            else:
                global_const.phis(j,k)=0.
            global_const.ps(j,k)=-RHOS*global_const.phis(j,k)-RHOS*u0*RE/2.*(2.*OM+u0/RE)*sn(k)**2
            u(j,k,:)=(global_const.ps(j,k)+global_const.ps0)*u0*cs(k)
    #
    #     Vorgabe einer Temperaturanomalie
    #
    for k in range(0,global_const.NK)
        for j in range(0,global_const.NJ)
            dist=sqrt((flam(j)-flamc)**2+(phi(k)-1.*phic)**2)
            if(dist.lt.distm):
                t(j,k,:)=(1.-dist/distm)
            else:
                t(j,k,:)=0.
            t(j,k,:)=(global_const.ps(j,k)+global_const.ps0)*t(j,k,:)