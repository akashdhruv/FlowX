def corrector(gridc, gridx, gridy, ivar, pvar, dt):


    u  = gridx.get_values(ivar)
    v  = gridy.get_values(ivar)    
    p  = gridc.get_values(pvar) 

    dx, dy = gridc.dx, gridc.dy

    u[1:-1,1:-1] = u[1:-1,1:-1] + dt*(p[2:-1,1:-1]-p[1:-2,1:-1])/dx
    v[1:-1,1:-1] = v[1:-1,1:-1] + dt*(p[1:-1,2:-1]-p[1:-1,1:-2])/dy

    gridx.fill_guard_cells(ivar)
    gridy.fill_guard_cells(ivar)

    return
