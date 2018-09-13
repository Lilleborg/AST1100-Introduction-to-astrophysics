def plot(pos,N,n):
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d.axes3d as p3
    import matplotlib.animation as animation
    """
    Plots the simulation of the particle Moving
    """
    def update_lines(num, dataLines,lines):
        for line, data in zip(lines,dataLines):
            line.set_data(data[0:2, num-1:num])
            line.set_3d_properties(data[2,num-1:num])
        return lines

    fig = plt.figure()
    ax = p3.Axes3D(fig)

    m = 100 #number of frames
    n = N  #number of particles
    N = n    #number of time steps

    data = pos
    #np.zeros([n,3,N])

    lines = [i for i in range(n)]
    for i in range(n):
        lines[i] = [ax.plot(data[i][0,0:1],
        data[i][1,0:1],data[i][2,0:1],'o')[0]]

    ax.set_xlim3d([0.0,L])
    ax.set_xlabel('X')

    ax.set_ylim3d([0.0,L])
    ax.set_ylabel('Y')

    ax.set_zlim3d([0.0,L])
    ax.set_zlabel('Z')

    ani = [i for i in range(n)]
    for i in range(n):
        ani[i] = animation.FuncAnimation(fig,update_lines,m,
        fargs=([data[i]],lines[i]),interval=50,blit=False)
    plt.show()
