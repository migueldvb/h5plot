#!/usr/bin/python2.5
"""
Calculate streamlines
"""

def stream(pxinit,pyinit,dirt):
# Read initial position and direction
    segment_nr = 5000
    iter_nr = 1000
    dt = 1
    drint = 0.01

    rx = xmin +(pxinit+0.5)*dx_fine
    ry = ymin +(pyinit+0.5)*dy_fine
    print "rx,ry = rx ry \n"

    xpos = (rx,rx)
    ypos = (ry,ry)

    for i in range(1,segment_nr):
        for j in range(1,iter_nr):
            get_velocity(rx,ry,vel)

            if (geometry eq "polar"):
                ry = ry + dirt*dt*vel[1]/rx
            else if (geometry eq "Cartesian"):
                ry = ry + dirt*dt*vel[1]
            rx = rx + dirt*dt*vel[0]
#             break if (( drint < sqrt((xpos[1]- rx)**2+(ypos[1]- ry)**2) && (geometry eq "Cartesian")))
#             break if (( drint < sqrt((xpos[1]- rx)**2+(xpos[1]*ypos[1]- rx*ry)**2) && (geometry eq "polar")))

        xpos[0] = xpos[1]
        xpos[1] = rx
        ypos[0] = ypos[1]
        ypos[1] = ry

        break if ((rx < xrange[0]) || (rx > xrange[1]))
        break if (((ry < yrange[0]) || (ry > yrange[1]-0.5*dy_fine)) && (geometry eq "Cartesian"))
#       last if (((ry < yrange[0]) || (ry > yrange[1]-0.5*dy_fine)) && (geometry eq "polar"))
#       last if ((0.5*drint) > sqrt((xpos[0]- rx)**2+(ypos[0]- ry)**2))
        break if (((0.5*drint) > sqrt((xpos[0]- rx)**2+(xpos[0]*ypos[0]- rx*ry)**2)) && (geometry eq "polar"))
        if ((ry < yrange[0]) && (geometry eq "polar")): ry = 2.*pi+ry
        if ((ry > yrange[1]) && (geometry eq "polar")): ry = ry-2.*pi

        if (! cart_mode): ypos = (ry,ry) 

        if (cart_mode):
            xcart[0] = xpos[0]*cos(ypos[0])
            ycart[0] = xpos[0]*sin(ypos[0])
            xcart[1] = xpos[1]*cos(ypos[1])
            ycart[1] = xpos[1]*sin(ypos[1])
            plt.plot(xcart,ycart)
        else:
            plt.plot(dist*xpos,dist*ypos)

def get_velocity(rxg,ryg,vel):
   ix = int((rxg - xrange[0] - 0.5*dx_fine)/dx_fine)# + 1
   qx = (rxg - xrange[0] - (ix+0.5)*dx_fine)/dx_fine

   iy = int(((ryg - yrange[0] - 0.5*dy_fine)/dy_fine))# + 1
   qy = (ryg - yrange[0] - (iy+0.5)*dy_fine)/dy_fine

   # Interpolate to calculate the velocity
   ix1 = ix+1
   iy1 = iy+1
   if ((iy >= ny) && (geometry eq "polar")): iy = iy%ny 
   if ((iy < 0) && (geometry eq "polar")): iy = iy+ny 
   if ((iy1 >= ny) && (geometry eq "polar")): iy1 = iy1%ny 
   if ((iy1 < 0) && (geometry eq "polar")): iy1 = iy1+ny 
   if ((ix >= nx) && (geometry eq "polar")): ix = ix%nx 
   if ((ix < 0) && (geometry eq "polar")): ix = ix+nx 
   if ((ix1 >= nx) && (geometry eq "polar")): ix1 = ix1%nx 
   if ((ix1 < 0) && (geometry eq "polar")): ix1 = ix1+nx 

#    print "nx,ny nx ny \n"
#    print "ix,iy ix iy \n"
#    print "ix1,iy1 ix1 iy1 \n"
   vel[0] = temp_velx(ix,iy) * (1-qx)*(1-qy) +
        temp_velx(ix1,iy) * qx*(1-qy) +
        temp_velx(ix1,iy1) * qx*qy +
        temp_velx(ix,iy1) * (1-qx)*qy

   vel[1] = temp_vely(ix,iy) * (1-qx)*(1-qy) +
        temp_vely(ix1,iy) * qx*(1-qy) +
        temp_vely(ix1,iy1) * qx*qy +
        temp_vely(ix,iy1) * (1-qx)*qy
#    print "vel vel \n"
