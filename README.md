# urdf_plot

A simple library that allows to visualize in Python a robot structure defined with URDF file. Example of usage: 

    rbt = Robot("my_robot.urdf")
    qs = {"joint1":0.7, "joint3":-0.3}
    rbt.plot(qs)