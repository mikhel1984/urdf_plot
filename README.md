# urdf_plot

A simple library that allows to visualize in Python a robot structure defined with URDF file. Example of usage: 

    rbt = Robot("my_robot.urdf")
    qs = {"joint1":0.7, "joint3":-0.3} # angles are optional
    rbt.plot(qs)

It can be used as a command line tool as well:

    ./urdfplot.py my_robot.urdf joint1 0.7 joint3 -0.3
    
